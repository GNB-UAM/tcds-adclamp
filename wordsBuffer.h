////////////////////////////////////////////////////////////////////////////////
/*!
 * @file wordsBuffer.h
 * @date 1-2013
 * @author Angel Lareo
 * 
 * @brief Definicion de la estructura WordsBuffer y sus funciones de operación
*/
////////////////////////////////////////////////////////////////////////////////

#ifndef __WORDS_BUFFER__
#define __WORDS_BUFFER__

#define MAX_WORDS_BUFF 60000

#define NOT_INITIATED -2
#define ERR -1
#define OK 0

#define MAX_BITS_WORD 10
#define MAX_WORDS 1024

typedef int WbElement;

/**
 * *struct BitsBuffer
 * @brief Definicion de la estructura BitsBuffer encargada de almacenar los bits
 * Funciona como un buffer circular de tamaño igual al de una palabra
 * El tamaño de palabra es definido por el usuario en GUI
 * 
 * @author Ángel Lareo
 * @date 1/2013
*/
typedef struct
{
    char bits[MAX_BITS_WORD];
    int wordLength;
    int numBits;
    char* init;
    char* insert;
} BitsBuffer;


/**
 * *struct WordsBuffer
 * @brief Definicion de la structura de WordsBuffer para almacenamiento de palabras
 * Array de palabras binarias almacenadas como enteros.
 * Consta de un BitsBuffer que almacena los bits de palabra temporalmente.
 * El máximo número de palabras viene definido por el usuario en GUI.
 * Si la ejecución del programa continua tras el almacenamiento del total de palabras,
 * las palabras iniciales se eliminarán secuencialmente.
 * 
 * @author Ángel Lareo
 * @date 12/2012
*/
typedef struct 
{
    WbElement words[MAX_WORDS_BUFF];
    BitsBuffer bb;
    int maxWords;
    int numWords;
    WbElement* init;
    WbElement* insert;
    WbElement* check;    
} WordsBuffer;

void wbCreateHistogram (WordsBuffer* wb, int* results, int numWords);
int wbInit(WordsBuffer* wb, int length, int maxWords);
int wbBitInsert(WordsBuffer* wb, char bit);  //Inserts bit on BitBuffer
int wbWordInsert (WordsBuffer* wb, int word);
int wbStoreWord (WordsBuffer *wb); //Transformation Bits 2 Double done here
int wbCheckWordMatch(WordsBuffer *wb, char *word);
int wbBits2Int(WordsBuffer *wb);
int Bits2Int(char* bb, int length); 

void bbAdvancePtr (BitsBuffer* bb, char** ptr);

#endif

