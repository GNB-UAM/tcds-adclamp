#include <stdio.h>
#include <stdlib.h>
#include "wordsBuffer.h"

#define NUM_OPTS 6
#define BITS_PER_WORD 5
#define MAX_WORDS 200

enum MenuOpts {EXIT, INIT, BIT_INSERT, WORD_INSERT, STORE_WORD, CHECK_WORD};
int calls[NUM_OPTS];
WordsBuffer wb;

int main(){
    int i, opt, ret=-1;   

    for (i=0; i<NUM_OPTS; i++) calls[i]=0;

    while (1){
        opt = menuOpts();
        switch (opt){
            case EXIT:
                calls[EXIT]++;
                return 0;
                break;
            case INIT:
                calls[INIT]++;
                ret = initTest();
                break;
            case BIT_INSERT:
                if (calls[INIT]==0) initTest();
                calls [BIT_INSERT]++;
                ret = bitInsertTest();
                break;
            case WORD_INSERT:
                if (calls[INIT]==0) initTest();
                calls [WORD_INSERT]++;
                wordInsertTest();
                break;
            case STORE_WORD:
                if (calls[INIT]==0) initTest();
                calls [STORE_WORD]++;
                storeWordTest();
                break;
            case CHECK_WORD:
                if (calls[INIT]==0){
					printf("Init\n"); initTest();
				}
                calls [CHECK_WORD]++;
                checkWordTest2();
                break;
            default: break;
        }

        printf ("Resultado del Test: %d %s\n", ret, ret==0?"OK":"ERROR"); 
    }
}

int menuOpts(){
    int out;
    printf ("\n----MENU----\n");
    printf ("%d) Salir\n", EXIT);
    printf ("%d) Iniciar wordsBuffer\n", INIT);
    printf ("%d) Insertar bits\n", BIT_INSERT);
    printf ("%d) Insertar word\n", WORD_INSERT);
    printf ("%d) Comprobar storeWord\n", STORE_WORD);
    printf ("%d) Comprobar wordMatch\n", CHECK_WORD);
    printf ("Introduzca opción: ");
    scanf("%d",&out);
    return out;
}

int initTest(){
    wbInit (&wb, BITS_PER_WORD, MAX_WORDS);
    return 0;
}

int bitInsertTest(){
    int i;
    char d;
    for (i=0; i<15; i++){
        d = i%2;
        wbBitInsert(&wb,d);
    }
    return 0;
}

int wordInsertTest(){
    int i;
    int d=1;
    WbElement *auxWordPtr;

    for (i=0; i<5; i++){
        wbWordInsert (&wb, d);
        d *= 2;
    }

    d = 1;
    auxWordPtr = wb.init;

    if (wb.numWords != 5) return -1;

    for (i=0; i<5; i++){
        if (auxWordPtr[i] != d) return -1;
        d *= 2;
    }

    return 0;
}

int storeWordTest(){
    int i;
    char d;

    for (i=0; i<1500; i++){
        d = i%2;
        wbBitInsert(&wb, d);
        if (i>=(BITS_PER_WORD-1)){
            wbStoreWord (&wb);
        }
    }

    if (wb.numWords != 1000) return -1;
    
    return 0;
}

int checkWordTest(){
    int i;
    char d, w[BITS_PER_WORD];
    for (i=0; i<BITS_PER_WORD; i++) w[i]=0;
    w[0]=1; w[BITS_PER_WORD-1]=1;

    /*for (i=0; i<1500; i++){
        wbBitInsert(&wb, d);
        if (i>=BITS_PER_WORD-1){
            wbStoreWord (&wb);
            if (1==wbCheckWordMatch(&wb, w)) return -1;
        }
    }*/

    for (i=0; i<BITS_PER_WORD; i++){
        wbBitInsert(&wb, w[i]);
    }
    wbStoreWord (&wb);

    printf("palabra: %d\n", wb.words[0]);

    return wbCheckWordMatch(&wb, w);
}

int checkWordTest2(){
    int i, detectedWord;
    char in = '0', w[BITS_PER_WORD];	
	
    for (i=0; i<BITS_PER_WORD; i++) w[i]=0;
    w[0]=1; w[BITS_PER_WORD-1]=1;
	
	fflush(stdin);
	getc(stdin);
	
	while (in!='-'){

		wbPrint(&wb, w);
		printf ("wordMatch: ");
			for (i=0; i<BITS_PER_WORD; i++) printf("%d", w[i]);
		printf ("\n");
		
		in = getc (stdin);
		in = atoi((char*) &in);
		
		getc (stdin);
		fflush(stdin);
		
		wbBitInsert(&wb,in);
		if (-1!=wbStoreWord(&wb)){
			detectedWord = wbCheckWordMatch(&wb, w);
			printf ("dw: %d\n", detectedWord);
		}
	}
	
    return wbCheckWordMatch(&wb, w);
}
