////////////////////////////////////////////////////////////////////////////////
/*!
 * @file words.c
 * @date 12-2012
 * @author Angel Lareo
 * 
 * @brief Archivo con las funciones del modelo words
 * 
 * Archivo con las funciones del modelo words. Este modelo detecta palabras en el
 * alfabeto binario estímulo/falta de estímulo y se encarga de disparar la estimulación
 * en tiempo real en función de dicha detección.
*/
////////////////////////////////////////////////////////////////////////////////

#include "words_shm.h"
#include "words.h"

#define __DEBUG__

//!< External Variables
NEURON_MODEL *neuronModel;
NEURON_MODEL *stimulator;

// To include a new variable that interacts with the GUI
// a new field should be added to the WORDS_SHM structure.
WORDS_SHM *wordsShm;  // Models variable
// ComediShm is used to interact with the DAQ board
DAQSHM *comediShm;

WordsBuffer wordsbuf;

extern TOPOLOGY *topology;

/*! init_neurons
   \description initializates the variables and parameters of a Neuron model comprises
   of several neurons
   \param NEURON_MODEL *neuronModel
   \return void
*/
void init_neurons (NEURON_MODEL *neuronModel)
{
    unsigned long i, pos_var=0, pos_param=0;

    DEBUG_LOG("Start init_neurons.\n");

    if ( !neuronModel )
        return;

    neuronModel->n_neurons = 1;
    for (i=0; i<neuronModel->n_neurons; i++)
    {
        //!< Set Variables
        neuronModel->variables[pos_var+X] = X_0;
        neuronModel->variables[pos_var+Y] = Y_0;
        neuronModel->variables[pos_var+Z] = Z_0;

        neuronModel->old_variables[pos_var+X] = neuronModel->variables[pos_var+X];
        neuronModel->old_variables[pos_var+Y] = neuronModel->variables[pos_var+Y];
        neuronModel->old_variables[pos_var+Z] = neuronModel->variables[pos_var+Z];

        //!< Pass V variable to FIFO
        set_nm_variable_fifo (neuronModel, pos_var, X);

        //!< Set Parameters
        neuronModel->params[pos_param+I_EXT] = I_0;
        neuronModel->params[pos_param+E] =  E_0;
        neuronModel->params[pos_param+K] =  K_0;

        //!< Pass I variable to FIFO
        //set_neuronModel_param_fifo (neuronModel, pos_param, I_EXT);

        //< Next Neuron
        pos_var += neuronModel->n_variables;
        pos_param += neuronModel->n_params;
    }

    if (neuronModel->n_neurons > 0)
        neuronModel->params[E] = 3.0;

    DEBUG_LOG("End init_neurons.\n");
    return;
}// End init_neurons


/*! asign_neurons_names
   \description Set neurons variables and parameters names
   \param NEURON_MODEL *neuronModel
   \return void
*/
void asign_neurons_names (NEURON_MODEL *neuronModel)
{
    if (!neuronModel)
        return;

    DEBUG_LOG ("Start asign_neurons_names.\n");

    //!< Neuron Model Name
    sprintf( neuronModel->neuron_model_name, "words" );

    //!< Set names to state variables
    set_nm_var_info (neuronModel, X, "X", MILIVOLTS, VOLTAGE, DOUBLE);
    set_nm_var_info (neuronModel, Y, "Y", NOUNITS, GATE_VARIABLE, DOUBLE);
    set_nm_var_info (neuronModel, Z, "Z", NOUNITS, GATE_VARIABLE, DOUBLE);

    //!<  Set names to parameters
    set_nm_param_info (neuronModel, I_EXT, "Iext", NANOAMPS, INTENSITY, DOUBLE);
    set_nm_param_info (neuronModel, E,  "E", NANOAMPS, INTENSITY, DOUBLE);
    set_nm_param_info (neuronModel, K,  "K", NOUNITS, SCALE, DOUBLE);

    DEBUG_LOG ("End asign_neurons_names.\n");

    return;
}// End of asign_neurons_names





/**
 * init_neuronModel
 * @description Initializes the neuron module
 * @param NEURON_MODEL *neuronModel
 * @return void
*/
NEURON_MODEL *init_nm ( unsigned long n_neurons )
{
    unsigned long i;

    DEBUG_LOG("\nStart init_neuronModel.\n");
    rt_printk ("Start words _neuronModel.\n");

    if (!( neuronModel = (NEURON_MODEL *) rtai_kmalloc(nam2num(WORDS_ADDR), sizeof(NEURON_MODEL))))
    {
        rt_printk("ERROR: Failure to initialize model shared memory.\n");
        return NULL;
    }

    if (!( wordsShm = (WORDS_SHM *) rtai_kmalloc(nam2num(WORDSSHM), sizeof(WORDS_SHM))))
    {
        rt_printk("ERROR: Failure to initialize words shared memory.\n");
        return NULL;
    }

    if (!( stimulator = (NEURON_MODEL *) rtai_kmalloc(nam2num("STISHM"), sizeof(NEURON_MODEL))))
    {
        rt_printk("ERROR: Failure to initialize model shared memory.\n");
        return NULL;
    }

    if (!(comediShm = (DAQSHM *) rtai_kmalloc (nam2num(COMEDI_ADDR), sizeof(DAQSHM))))
    {
        rt_printk("ERROR: Failure to initialize shared memory.\n");
        return NULL;
    }


    neuronModel->n_neurons=n_neurons;
    neuronModel->n_variables=DIM;
    neuronModel->n_params=DIM_PARAM;
    neuronModel->n_variables_fifo=0;
    neuronModel->n_params_fifo=0;

//    rt_printk ("Assigning names.\n");

    //!< Assign names to neurons
    asign_neurons_names ( neuronModel );

//    rt_printk ("Initializing neurons.\n");
    //!< Initialize neurons
    init_neurons (neuronModel);

    for (i=0; i < neuronModel->n_neurons; i++)
        neuronModel->n_synapses[i]=0;

    neuronModel->integrate_neuron_model = nm_func2;
    neuronModel->trigger_function = 0;
    
    wordsShm->channel_in = 0; //Select channel input

    DEBUG_LOG("End init_neuronModel.\n");

    return neuronModel;
}// End init_neuronModel

/**
 * free_neuronModel
 * @description Initializes the neuron module
 * @param NEURON_MODEL *neuronModel
 * @return void
*/
void free_nm (void)
{
    rtai_kfree(nam2num(WORDS_ADDR));
    rtai_kfree(nam2num(WORDSSHM));
    rtai_kfree(nam2num("STISHM"));
    rtai_kfree(nam2num(COMEDI_ADDR));
}

/** int_function
 * @description This function is called by the integration method
 * @return void
*/
void int_function(double time, double *vars, double * dvars, double *p)
{
    // It is compulsary to declare, but not to fill in !!!!!!!!!!!
}

/** int_function
 * @description This function is called by the integration method
 * @return void
*/
void nm_func2(void)// This function is called in every step of the dynamic clamp
{
	neuronModel->variables[X] = model(comediShm->buffer_in[wordsShm->channel_in], topology->t);
	
    // model function is the one to be modified
    // time is in ms: topology->t
}

double model(double currentV, double time){
	if (wordsShm->activated[0]) return detectorModel(currentV, time);
	else if (wordsShm->activated[1]) return histogramModel(currentV, time);
	else if (wordsShm->activated[2]) return aleatModel(currentV, time);
}

double histogramModel(double currentV, double time){
	static double lastTime, initTime;
    static short numCalls = FIRST_CALL;
    static char out = BIT_NOT_DETECTED_OUT;
    static double bitTime;
    static WordsBuffer *wb;
	static int n=0;
	
    if (numCalls==FIRST_CALL){
		wb = &wordsbuf;
        wbInit(wb,wordsShm->wordSize, wordsShm->numWords);
		initTime = time;
        numCalls = NEW_TIME_WINDOW;
        n++;
        lastTime = time;
        bitTime = wordsShm->bitTime;
        numCalls++;
    }
    
	if (time-initTime < (wordsShm->windowTime*1000)){
		if ((detect_spike(currentV, time))&&(out != BIT_DETECTED_OUT)){ //DETECT SPIKE
			out = BIT_DETECTED_OUT;    
			wbBitInsert(wb,1); //Store 1 for this window time	
			wbStoreWord(wb);
		}

		if (time >= initTime+n*wordsShm->bitTime){ //Window time ended     
			if (out == BIT_NOT_DETECTED_OUT){
				wbBitInsert(wb,0); //Store 0 for this window time
				wbStoreWord(wb);
			}
			out = BIT_NOT_DETECTED_OUT;
            n++;
            lastTime = time;
            bitTime = wordsShm->bitTime;
            numCalls++;
			return END_WINDOW_OUT;
		}
	}
	else {
		initTime = time;
		wbCreateHistogram(wb, wordsShm->hist, wordsShm->numWords);
		wordsShm->histSem++;
		n=0;
		return END_HISTOGRAM_TIME;	
	}
	
	return out;
}

/*double detectorModel(double currentV, double time)
{
    static double lastTime;
    static short numCalls = FIRST_CALL;
    static char out = BIT_NOT_DETECTED_OUT, detectedWord = FALSE;
    static double bitTime;
	static WordsBuffer *wb;
	
    if (numCalls==FIRST_CALL){
		wb = &wordsbuf;
        wbInit(wb,wordsShm->wordSize, wordsShm->numWords);
        lastTime = time;
        bitTime = wordsShm->bitTime;
        numCalls++;
    }

    if ((detect_spike(currentV, time))&&(out != 1)){ //DETECT SPIKE
        out = BIT_DETECTED_OUT;    
        wbBitInsert(wb,out); //Store 1 for this window time
        if (-1!=wbStoreWord(wb)){
            detectedWord = wbCheckWordMatch(wb, wordsShm->word);
		}
    }

    if (bitTime < time-lastTime){ //Window time ended     
        if (out == BIT_NOT_DETECTED_OUT){
            wbBitInsert(wb,out); //Store 0 for this window time
            if (-1!=wbStoreWord(wb))
                detectedWord = wbCheckWordMatch(wb, wordsShm->word);
        }
        
        out = BIT_NOT_DETECTED_OUT;
        lastTime = time;
        bitTime = wordsShm->bitTime;
        numCalls++;
        
        if (detectedWord){
            if (stimulator->trigger_function){
				#ifdef __DEBUG__
				rt_printk ("\nStimulator trigger_function");
				#endif
				stimulator->trigger_function();
			}
			detectedWord = FALSE;
            return WORD_DETECTED_OUT;
        }

        detectedWord = FALSE;
        return END_WINDOW_OUT;
    }

    return out;
}*/

double detectorModel(double currentV, double time)
{
    static double lastTime;
    static short numCalls = FIRST_CALL;
    static char out = BIT_NOT_DETECTED_OUT, detectedWord = FALSE, fixedDelay=TRUE;
    static double bitTime, wordDetectTime=100000, stimTime= -400;
	static WordsBuffer *wb;

    
    if (numCalls==FIRST_CALL){
        if (wordsShm->minDelay<0) // Negative Delay = Not Fixed Delay
            fixedDelay=FALSE;
        //wordsShm->refractTime = 200;
		wb = &wordsbuf;
        wbInit(wb,wordsShm->wordSize, wordsShm->numWords);
        numCalls = NEW_TIME_WINDOW;
        lastTime = time;
        bitTime = wordsShm->bitTime;
        numCalls++;
    }

    if ((detect_spike(currentV, time))&&(out != 1)){ //DETECT SPIKE
        out = BIT_DETECTED_OUT;    
        wbBitInsert(wb,out); //Store 1 for this window time
        if (-1!=wbStoreWord(wb)){
            detectedWord = wbCheckWordMatch(wb, wordsShm->word);
            if (detectedWord) wordDetectTime=time;
		}
    }

    if (detectedWord && fixedDelay && time >= wordDetectTime + wordsShm->minDelay){ //precise delay
            if (stimulator->trigger_function){ 
            //if (stimulator->trigger_function && ((time - stimTime) > 400)){ //refract time
				#ifdef __DEBUG__
				rt_printk ("\nStimulator trigger_function");
				#endif
				stimulator->trigger_function();
                stimTime = time;
			}
			detectedWord = FALSE;
            return WORD_DETECTED_OUT;
    }

    if (bitTime < time-lastTime){ //Window time ended     
        if (out == BIT_NOT_DETECTED_OUT){
            wbBitInsert(wb,out); //Store 0 for this window time
            if (-1!=wbStoreWord(wb)){
                detectedWord = wbCheckWordMatch(wb, wordsShm->word);
                if (detectedWord) wordDetectTime=time;
            }
        }
        
        if (detectedWord && !fixedDelay){ //precise delay
            if (stimulator->trigger_function){ 
            //if (stimulator->trigger_function && ((time - stimTime) > 400)){ //refract time
				#ifdef __DEBUG__
				rt_printk ("\nStimulator trigger_function");
				#endif
				stimulator->trigger_function();
                stimTime = time;
			}
			detectedWord = FALSE;
            return WORD_DETECTED_OUT;
        }

        out = BIT_NOT_DETECTED_OUT;
        lastTime = time;
        numCalls++;

        return END_WINDOW_OUT;
    }

    return out;
}


/*double aleatModel (double voltage, double time){
    static double lastTime;
    static short numCalls = NEW_TIME_WINDOW;
	static unsigned int jitter, lanzado = 0;
	unsigned int i;
	
    if (numCalls==NEW_TIME_WINDOW){ //New time window
		lanzado = 0;
		get_random_bytes(&i, sizeof(unsigned int));
		jitter = i % (unsigned int)(wordsShm->windowTime*1000);
        rt_printk ("\nJitter: %u", jitter);
        lastTime = time;
        numCalls++;
    }
    
	if ((jitter < (unsigned int)(time-lastTime)) && (lanzado == 0)){
		if (stimulator->trigger_function){
			#ifdef __DEBUG__
			rt_printk ("\nInfo stimulator trigger_function");
			#endif
			stimulator->trigger_function();
		}
		lanzado = 1;
	}
	
	if ((wordsShm->windowTime*1000) < time-lastTime){
		numCalls = NEW_TIME_WINDOW;
	}
	
	return 0;
}*/

/*double aleatModel (double voltage, double time){
    static double lastTime, spkTime;
    static short numCalls = NEW_TIME_WINDOW;
	static unsigned int jitter, jitter2, lanzado = 0, spkDetected = 0;
	unsigned int i;
    
    if (stimulator->trigger_function){
    	spkDetected = detect_spike(voltage, time);
        if (spkDetected==1) spkTime = time;

        if (numCalls==NEW_TIME_WINDOW){ //New period
		    lanzado = 0;
            spkDetected = 0;
		    get_random_bytes(&i, sizeof(unsigned int));
		    jitter = i % (unsigned int)(wordsShm->randomPeriod*1000 - wordsShm->windowTime);
            lastTime = time;
            numCalls++;
            get_random_bytes(&i, sizeof(unsigned int));
            jitter2 = (wordsShm->wordSize)*(wordsShm->windowTime) + (i % (unsigned int)(wordsShm->windowTime));
        }

        if ((lanzado==0) && (jitter < (unsigned int)(time-lastTime))){ //jitter passed
            if ((lanzado == 0) && (jitter2 < (unsigned int)(time - spkTime)) &&
                          (((wordsShm->wordSize)*(wordsShm->windowTime) + wordsShm->windowTime) > (time - spkTime))){
		        #ifdef __DEBUG__
		        rt_printk ("Info stimulator trigger_function %u < %u\n", jitter2, (unsigned int)(time - spkTime));
		        #endif
		        stimulator->trigger_function();
                lanzado=1;
            }
	    }
	
        if ((wordsShm->randomPeriod*1000) < time-lastTime){ //randomPeriod passed
	        numCalls = NEW_TIME_WINDOW;
      	}
    }
	
	return 0;
}*/

double aleatModel (double voltage, double time){
    static double lastTime, spkTime;
    static short numCalls = NEW_TIME_WINDOW;
	static unsigned int jitter, jitter2, lanzado = 0, spkDetected = 0;
	unsigned int i;
    
    if (stimulator->trigger_function){
    	spkDetected = detect_spike(voltage, time);
        if (spkDetected==1) spkTime = time;

        if (numCalls==NEW_TIME_WINDOW){ //New period
		    lanzado = 0;
            spkDetected = 0;
		    get_random_bytes(&i, sizeof(unsigned int));
		    jitter = i % (unsigned int)(wordsShm->randomPeriod*1000 - wordsShm->maxDelay);
            lastTime = time;
            numCalls++;
            get_random_bytes(&i, sizeof(unsigned int));
            if (wordsShm->minDelay == wordsShm->maxDelay)
                jitter2 = wordsShm->minDelay;
            else
                jitter2 = (wordsShm->minDelay) + (i % (unsigned int)(wordsShm->maxDelay - wordsShm->minDelay));
        }

        if ((lanzado==0) && (jitter < (unsigned int)(time-lastTime))){ //jitter passed
            if ((lanzado == 0) && (jitter2 < (unsigned int)(time - spkTime))){
		        #ifdef __DEBUG__
		        rt_printk ("Info stimulator trigger_function %u < %u\n", jitter2, (unsigned int)(time - spkTime));
		        #endif
		        stimulator->trigger_function();
                lanzado=1;
            }
	    }
	
        if ((wordsShm->randomPeriod*1000) < time-lastTime){ //randomPeriod passed
	        numCalls = NEW_TIME_WINDOW;
      	}
    }
	
	return 0;
}

int detect_spike (double voltage, double time)
{
    static double lastVolt2=0.0, lastVolt1=0.0;
    double volThreshold = wordsShm->threshold;
    
    int spike=0;
    
    if ((lastVolt2<=lastVolt1) &&
        (lastVolt1>=voltage) &&
        (lastVolt1>volThreshold))
    {
        spike=1;
    }

    lastVolt2 = lastVolt1;
    lastVolt1 = voltage;

    return spike;
}
