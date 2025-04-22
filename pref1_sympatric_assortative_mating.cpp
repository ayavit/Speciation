// Diploid population. BDM with identified mutations. Dominant effect (incompatibility due to one chromosom results in inviability).

#include <stdio.h>
#include <stdlib.h>
#include "math.h"
// #using <mscorlib.dll>
#include <search.h>
#include <time.h>
// using namespace System;


// probability that a couple of mutations will be incompatible
// This is the same probability for two derived mutations and for one derived and one ancestral.
double CONT_PROB=1./10;
clock_t t1, t2, t3, t4, tot1=0, tot2=0, tot3=0,tot4=0,tot11=0;
#define MASK 0XFFFFUL
unsigned long saf = CONT_PROB*(MASK+1);

// type of structure that contains a pointer to an array of the individual's mutations, a pointer to a boolean array that has TRUE where
// the individual has the mutation twice, and the number of mutations the individual has
typedef struct {
	unsigned long* list;  // pointer to an array of mutations that will be allocated for the individual
	unsigned long* list_hash;  // pointer to an array of hash values of mutations that will be allocated for the individual
	bool* two_copies;  // pointer to a boolean array that has TRUE where the individual has the mutation twice
	unsigned long n;  // length of the array. The number of different mutations the individual has.
} mut_st;


char* Itoa(int value, char* str, int radix) {
    static char dig[] =
        "0123456789"
        "abcdefghijklmnopqrstuvwxyz";
    int n = 0, neg = 0;
    unsigned int v;
    char* p, *q;
    char c;

    if (radix == 10 && value < 0) {
        value = -value;
        neg = 1;
    }
    v = value;
    do {
        str[n++] = dig[v%radix];
        v /= radix;
    } while (v);
    if (neg)
        str[n++] = '-';
    str[n] = '\0';

    for (p = str, q = p + (n-1); p < q; ++p, --q)
        c = *p, *p = *q, *q = c;
    return str;
}


int comp (const void *e1, const void *e2) 
{
	return *((int *) e1) - *((int *)e2);
}


long hash(int i1)
{
	return 32347*((unsigned long)i1+12343)+ 42347*(unsigned long)(i1+12433);
}

// function that returns false if the couple of mutations is compatible, and true if they are incompatible
bool contradict (unsigned long i1, unsigned long i2)
{
//	t1 = clock();
	if (i1==i2) return false;

	// very big number that is a function of the two mutation numbers
	unsigned long hash_val = i1 + i2;
//		32347*((unsigned long)i1+12343) + 42347*(unsigned long)(i2+12433)
		//+ 32347*((unsigned long)i2+12343) + 42347*(unsigned long)(i1+12433);

//	hash_val = hash_val % MASK + 1;
	hash_val = hash_val & MASK;
//	tot11 += clock()-t1;

	if (hash_val < saf) {/*tot1 += clock()-t1;*/ 
		return true;}
//	tot1 += clock()-t1;
	return false;
}

// function that returns false if the couple of mutations is compatible, and true if they are incompatible
bool contradict_try (unsigned long i1, unsigned long i2)
{
//	t1 = clock();
	if (i1==i2) return false;

	// very big number that is a function of the two mutation numbers
	unsigned long hash_val = 32347*((unsigned long)i1+12343) + 42347*(unsigned long)(i2+12433)
		+ 32347*((unsigned long)i2+12343) + 42347*(unsigned long)(i1+12433);

	hash_val = hash_val % RAND_MAX;
	if (hash_val < CONT_PROB * RAND_MAX) {/*tot1 += clock()-t1;*/ return true;}
//	tot1 += clock()-t1;
	return false;
}


// function that returns false if the couple of mutations is compatible, and true if they are incompatible
bool contradict_old (unsigned int i1, unsigned int i2)
{
	t1 = clock();
	if (i1==i2) return false;

	// very big number that is a function of the two mutation numbers
	unsigned long hash_val = 32347*((unsigned long)i1+12343) + 42347*(unsigned long)(i2+12433)
		+ 32347*((unsigned long)i2+12343) + 42347*(unsigned long)(i1+12433);

	hash_val = hash_val % RAND_MAX;
	if (hash_val < CONT_PROB * RAND_MAX) {tot1 += clock()-t1; return true;}
	tot1 += clock()-t1;
	return false;
}



// function that receives pointers to two mutation structures and returns 1 if the mating of the two was successful (the offspring is viable)
// and 0 if not. In addition, it returns the (first) mutations that caused inviability.
int is_success (mut_st *mp1_point, mut_st *mp2_point, unsigned int *cause1, unsigned int *cause2)
{

    
	// pointers to the first two causative mutations
	(*cause1) = 0;
	(*cause2) = 0;

	// the mutation structures: mp1 and mp2 
	mut_st mp1 = *mp1_point, mp2 = *mp2_point;

	// If the parents have no mutations, the mating is successful.
	if (mp1.n==0 && mp2.n==0) return 1;

    bool *p1_m1=NULL, *p1_m2=NULL;
//	t2 = clock();
	
    // allocate a boolean for every mutation
	p1_m1 = (bool *) malloc(mp1.n*sizeof(bool));  // indicator for getting a mutation from parent 1
	p1_m2 = (bool *) malloc(mp2.n*sizeof(bool));  // indicator for getting a mutation from parent 2

	// Each mutation from parent 2 is passed to the offspring with probability 0.5 or 1.
	for (int m2 = 0; m2<mp2.n ; m2++){
		if (mp2.two_copies[m2])
			p1_m2[m2] = true;
		else
			p1_m2[m2] = rand() >= RAND_MAX/2;
	}

	// Each mutation from parent 1 is passed to the offspring with probability 0.5 or 1.
	for (int m1 = 0; m1<mp1.n ; m1++){
		if (mp1.two_copies[m1])
			p1_m1[m1] = true;
		else
			p1_m1[m1] = rand() >= RAND_MAX/2;
	}


	// new offspring
    mut_st new_record;
	new_record.n = 0;
	new_record.list = NULL;
	new_record.list_hash = NULL;
	new_record.two_copies = NULL;

	// count the (maximal number of) new offspring's different mutations
	for (int m2 = 0; m2<mp2.n ; m2++) if (p1_m2[m2]) new_record.n++;
	for (int m1 = 0; m1<mp1.n ; m1++) if (p1_m1[m1]) new_record.n++;

	// create list of offspring mutations
	// if (new_record.n>0) 
	new_record.list = (unsigned long *) calloc (new_record.n+1,sizeof(long));
	new_record.list_hash = (unsigned long *) calloc (new_record.n+1,sizeof(long));
	new_record.two_copies = (bool *) calloc (new_record.n+1,sizeof(bool));
	
	// record the offspring's mutations
	
	int temp = 0;   // index in offspring's mutation list
	int i1 = 0;   // index in parent 1's mutation list
	int i2 = 0;   // index in parent 2's mutation list

	while (true){

		// stop if we reached the end of a list
		if (i1 == mp1.n && i2 == mp2.n) break;

		// if the mutation number in place i1 in first parent's list is bigger than the mutation number in place i2 in second
		// parent's list, record the mutation in index i2 (if passed to the offspring) and move index i2 forward.
		if ((i1 == mp1.n)||(i2 < mp2.n && mp1.list[i1]>mp2.list[i2])){

		    if (p1_m2[i2]){
		        unsigned int mm2 = mp2.list[i2];
		        new_record.list[temp] = mm2;
			    new_record.list_hash[temp] = mp2.list_hash[i2];
			    new_record.two_copies[temp] = false;
			    temp++;
			}

			i2++;
	    }

		else{

	        // if the mutation number in place i1 in first parent's list is smaller than the mutation number in place i2 in second
	        // parent's list, record the mutation in index i1 (if passed to the offspring) and move the index i1 forward.
	        if ((i2 == mp2.n)||(i1 < mp1.n && mp1.list[i1]<mp2.list[i2]) ){

	            if (p1_m1[i1]){
		            unsigned int mm1 = mp1.list[i1];
		            new_record.list[temp] = mm1;
		            new_record.list_hash[temp] = mp1.list_hash[i1];
		            new_record.two_copies[temp] = false;
		            temp++;
	            }

				i1++;
	        }

			else{

	            // If we found the same mutation in both parents, record the mutation and the number of copies passes to the offspring (if passed) and move indices forward.
	            // Also, subtract 1 from the offspring's count different mutations. 
	            if (i1 < mp1.n && i2 < mp2.n &&mp1.list[i1]==mp2.list[i2]){
	                if (p1_m1[i1] && p1_m2[i2]){
		                unsigned int mm1 = mp1.list[i1];
		                new_record.list[temp] = mm1;
		                new_record.list_hash[temp] = mp1.list_hash[i1];
		                new_record.two_copies[temp] = true;
						new_record.n--;
		                temp++;
	                }
	                else{
	                    if (p1_m1[i1]){
	                        unsigned int mm1 = mp1.list[i1];
		                    new_record.list[temp] = mm1;
		                    new_record.list_hash[temp] = mp1.list_hash[i1];
		                    new_record.two_copies[temp] = false;
		                    temp++;
		                }
		                if (p1_m2[i2]){
		                    unsigned int mm2 = mp2.list[i2];
		                    new_record.list[temp] = mm2;
		                    new_record.list_hash[temp] = mp2.list_hash[i2];
		                    new_record.two_copies[temp] = false;
		                    temp++;
		                }
		            }

		            i1++;
					i2++;
		        }
			}
		}
	}
	
	// check for incompatibility
	for (int m1 = 0; m1<new_record.n ; m1++){
		if (new_record.two_copies[m1] == false){
            for (int m2 = (m1+1); m2<new_record.n ; m2++){
			    if(new_record.two_copies[m2] == false){
					if (contradict(new_record.list_hash[m1],new_record.list_hash[m2])){
				
                        // save the mutations that caused inviability
				        (*cause1) = new_record.list[m1];
				        (*cause2) = new_record.list[m2];

				        goto RET0;

			        }
		        }
		    }
		}
	}


	free(p1_m1);
	free(p1_m2);
	free(new_record.list);
	free(new_record.list_hash);
	free(new_record.two_copies);

	return 1;


    RET0 :

	   free(p1_m1);
	   free(p1_m2);
	   free(new_record.list);
	   free(new_record.list_hash);
	   free(new_record.two_copies);
	
	   return 0;

}





int main(int argc, char* argv[])
{
		
	int NCH = 1000;   //size of population 
	int NGEN = 150000;   // number of generations to run
	double MUT_PROB = 0.1;   // probability of new mutation in a new individual
    double theta = 0;   // exponential function parameter. If theta=0 then the choice of parents does not depend on the distance between them.
	
	unsigned int seed_for_rand = 2;   // seed for randomization
    
	unsigned int inv_cause1 = 0, inv_cause2 = 0;   // the (first) mutations that cause inviability of a couple's offspring
	char mut_string1 [33];   // string to convert the first causative mutation number into
	char mut_string2 [33];   // string to convert the second causative mutation number into
	char success_string [33];   // string to convert the number of successes into
	char fname_all_mutations[200];   // string to write the name of the file to write all the individuals' mutations into
	char individual_number_string[33];   // string to contain an individual's number
	unsigned int number_of_mutations;   // number of mutations an individual has
	char this_mutation_string[33];   // a mutation an individual has, converted to string
    char fname[200];   // names of files
	char ffname[200];   // names of files
	char fffname[200];   // names of files
	char ffffname[200];   // names of files
	FILE *f;   // pointer to file
    FILE *ff;   // pointer to file 
	FILE *fff;   // pointer to file
	FILE *ffff;   // pointer to file

	// the user may deliver parameters to the program
	if (argc >= 2) NCH = atoi ((char *)(argv[1]));
	if (argc >= 3) CONT_PROB = atof ((char *)(argv[2])) ;
	if (argc >= 4) NGEN = atoi ((char *)(argv[3]));
	if (argc >= 5) MUT_PROB = atof ((char *)(argv[4]));
	if (argc >= 6) theta = atof ((char *)(argv[5]));
	if (argc >= 7) seed_for_rand = atoi ((char *)(argv[6]));

	// print the parameters
	printf ("NCH %d, CONT_PROB %f, NGEN %d, MUT_PROB %f, theta %f, seed_for_rand %d \n",NCH, CONT_PROB, NGEN, MUT_PROB, theta, seed_for_rand);

	// distances between individuals to be checked for compatibility in mating
	int try_dists[5] = {1,10,50,100,500};

	// list of all mutation structures of the individuals
	mut_st *mut_list;
	mut_st *new_mut_list;

	// allocate an array of NCH mutation structures 
	mut_list = (mut_st *) calloc (NCH, sizeof(mut_st));
	new_mut_list = (mut_st *) calloc (NCH, sizeof(mut_st));

	// set the seed for the rest of the rand() calls in the program
	srand(seed_for_rand);

	int total_mut = 0;   // total mutations counter
	mut_st mp1, mp2, new_record;   // mutation structures of two parents and offspring
	bool alive;   // is the offspring alive?
	bool *p1_m1 = NULL, *p1_m2 = NULL;
	bool *dup1 = NULL, *dup2 = NULL;
//	clock_t t1, t1_tot = 0, t2, t2_tot = 0, t3, t3_tot = 0, t4, t4_tot = 0;
	int sgn1, sgn2;
	double rn1, rn2;
	int tmp1, tmp2;
	int p1, p2;
	int success;
	int failures;
	int gen, i, t, r;   // loop indices

    


	// main loop: going over all the generations
	for (gen=0; gen<NGEN; gen++){

		failures = 0;   // counting failures (number of inviable individuals produced) in the current generation


		// create NCH (size of population) new individuals
		for (i=0; i<NCH; i++){

			failures--;   // If we got here, then the last failure was actually a success.  
			
			// pick parent 1
			sgn1 = (rand() > RAND_MAX/2) ? 1 : -1;	

            // if theta=0, choose parent 1 randomly
			if (theta==0){
					tmp1 = (((unsigned long)rand()*0x10000) + rand()) % NCH;
					//tmp2 = (((unsigned long)rand()*0x10000) + rand()) % NCH;
			}
								
			// if theta>0, choose parent 1 using exponential decay function
			else
			    do {
				    rn1 = ((double) rand())/RAND_MAX+0.00000001;
				    tmp1 = (int) (-log(rn1)/theta);
			    } while(tmp1>=NCH/2);

			// parent 1 chosen
			p1 = (i + sgn1 * tmp1 + NCH) % NCH ;
			//p2 = (i + sgn2 * tmp2 + NCH) % NCH ;

			// mutation structure of parent 1
			mp1 = mut_list[p1];
			//mp2 = mut_list[p2];

			bool *p1_m1=NULL;

            // allocate a boolean for every mutation
	        p1_m1 = (bool *) malloc(mp1.n*sizeof(bool));  // indicator for getting a mutation from parent 1
	        //p1_m2 = (bool *) malloc(mp2.n*sizeof(bool));  // indicator for getting a mutation from parent 2


			// produce a viable offspring
			do {

                // initialize new offspring
                new_record.n = 0;
	            new_record.list = NULL;
	            new_record.list_hash = NULL;
			    new_record.two_copies = NULL;

				
				failures++;   // increase the failures' counter anyway
				alive = true;   // unless we find an incompatibility, assume the offspring is viable

                
				// pick parent 2

				//sgn1 = (rand() > RAND_MAX/2) ? 1 : -1;
				sgn2 = (rand() > RAND_MAX/2) ? 1 : -1; 

				// if theta=0, choose parent 2 randomly
				if (theta==0){
					//tmp1 = (((unsigned long)rand()*0x10000) + rand()) % NCH;
					tmp2 = (((unsigned long)rand()*0x10000) + rand()) % NCH;
				}
								
				// if theta>0, choose parent 2 using exponential decay function
				else
				    do {
					    rn2 = ((double) rand())/RAND_MAX+0.00000001;
					    tmp2 = (int) (-log(rn2)/theta);
				    } while(tmp2>=NCH/2);

				// parent 2 chosen
				//p1 = (i + sgn1 * tmp1 + NCH) % NCH ;
				p2 = (i + sgn2 * tmp2 + NCH) % NCH ;


				// mutation structure of parent 2
				//mp1 = mut_list[p1];
				mp2 = mut_list[p2];

				bool *p1_m2=NULL;

                // allocate a boolean for every mutation
	            //p1_m1 = (bool *) malloc(mp1.n*sizeof(bool));  // indicator for getting a mutation from parent 1
	            p1_m2 = (bool *) malloc(mp2.n*sizeof(bool));  // indicator for getting a mutation from parent 2

	            // Each mutation from parent 2 is passed to the offspring with probability 0.5 or 1.
	            for (int m2 = 0; m2<mp2.n ; m2++){
		            if (mp2.two_copies[m2])
			            p1_m2[m2] = true;
		            else
			            p1_m2[m2] = rand() >= RAND_MAX/2;
	            } 

	            // Each mutation from parent 1 is passed to the offspring with probability 0.5 or 1.
	            for (int m1 = 0; m1<mp1.n ; m1++){
		            if (mp1.two_copies[m1])
			            p1_m1[m1] = true;
		            else
			            p1_m1[m1] = rand() >= RAND_MAX/2;
	            }



	            // count the (maximal number of) new offspring's different mutations
	            for (int m2 = 0; m2<mp2.n ; m2++) if (p1_m2[m2]) new_record.n++;
	            for (int m1 = 0; m1<mp1.n ; m1++) if (p1_m1[m1]) new_record.n++;

	            // create list of offspring mutations
	            // if (new_record.n>0) 
	            new_record.list = (unsigned long *) calloc (new_record.n+1,sizeof(long));
	            new_record.list_hash = (unsigned long *) calloc (new_record.n+1,sizeof(long));
				new_record.two_copies = (bool *) calloc (new_record.n+1,sizeof(bool));

	
	            // record the offspring's mutations
	
	            int temp = 0;   // index in offspring's mutation list
	            int i1 = 0;   // index in parent 1's mutation list
	            int i2 = 0;   // index in parent 2's mutation list

	            while (true){

		            // stop if we reached the end of a list
		            if (i1 == mp1.n && i2 == mp2.n) break;

		            // if the mutation number in place i1 in first parent's list is bigger than the mutation number in place i2 in second
		            // parent's list, record the mutation in index i2 (if passed to the offspring) and move index i2 forward.
		            if ((i1 == mp1.n)||(i2 < mp2.n && mp1.list[i1]>mp2.list[i2])){

			            if (p1_m2[i2]){
					       unsigned int mm2 = mp2.list[i2];
					       new_record.list[temp] = mm2;
					       new_record.list_hash[temp] = mp2.list_hash[i2];
					       new_record.two_copies[temp] = false;
					       temp++;
			            }

			            i2++;
			        }

					else{

		                // if the mutation number in place i1 in first parent's list is smaller than the mutation number in place i2 in second
		                // parent's list, record the mutation in index i1 (if passed to the offspring) and move the index i1 forward.
		                if ((i2 == mp2.n)||(i1 < mp1.n && mp1.list[i1]<mp2.list[i2]) ){

			                if (p1_m1[i1]){
					            unsigned int mm1 = mp1.list[i1];
					            new_record.list[temp] = mm1;
					            new_record.list_hash[temp] = mp1.list_hash[i1];
					            new_record.two_copies[temp] = false;
					            temp++;
			                }

							i1++;
			            }

						else{

		                    // If we found the same mutation in both parents, record the mutation and the number of copies passed to the offspring (if passed) and move indices forward.
		                    // Also, subtract 1 from the offspring's count ofdifferent mutations. 
		                    if (i1 < mp1.n && i2 < mp2.n &&mp1.list[i1]==mp2.list[i2]){
			                    if (p1_m1[i1] && p1_m2[i2]){
					                unsigned int mm1 = mp1.list[i1];
					                new_record.list[temp] = mm1;
					                new_record.list_hash[temp] = mp1.list_hash[i1];
					                new_record.two_copies[temp] = true;
									new_record.n--;
					                temp++;
			                    }
			                    else{
				                    if (p1_m1[i1]){
					                    unsigned int mm1 = mp1.list[i1];
					                    new_record.list[temp] = mm1;
					                    new_record.list_hash[temp] = mp1.list_hash[i1];
					                    new_record.two_copies[temp] = false;
					                    temp++;
				                    }
				                    if (p1_m2[i2]){
					                    unsigned int mm2 = mp2.list[i2];
					                    new_record.list[temp] = mm2;
					                    new_record.list_hash[temp] = mp2.list_hash[i2];
					                    new_record.two_copies[temp] = false;
					                    temp++;
			                        }
			                    }

			                    i1++;i2++;
		                    }
						}
					}
	            }


	            // check for incompatibility
	            for (int m1 = 0; m1<new_record.n ; m1++)
		            if (new_record.two_copies[m1] == false)
                        for (int m2 = (m1+1); m2<new_record.n ; m2++)
			                if(new_record.two_copies[m2] == false)
					            if (contradict(new_record.list_hash[m1],new_record.list_hash[m2]))
						            alive=false;
	
	            
	            free(p1_m2);

				if (!alive){
                    free (new_record.list);
			        free (new_record.list_hash);
			        free (new_record.two_copies);
				}


            // continue until a viable offspring is produced
			} while (!alive);
			
			free(p1_m1);
			

			// replace individual i in the population with the new offspring
//			if (gen > 0){
//				free (mut_list[i].list);
//				free (mut_list[i].list_hash);
//				free (mut_list[i].two_copies);
//			}
			new_mut_list[i] = new_record;

			//free (new_record.list);
			//free (new_record.list_hash);
			//free (new_record.two_copies);

			//qsort (mut_list[i].list,  mut_list[i].n, sizeof(int), comp);

		}   // end of loop for creating NCH new offsprings
		


		for (i = 0 ; i<NCH ; i++){
			if (gen > 0){
				free (mut_list[i].list);
				free (mut_list[i].list_hash);
				free (mut_list[i].two_copies);
			}
			mut_list[i] = new_mut_list[i];
		}
		
		
		// development of new mutations in the population in this generation
		for (i=0; i<NCH; i++){

			// Each individual has MUT_PROB probability to develop a new mutation. 
			if (rand()< MUT_PROB * RAND_MAX){
				total_mut++;
				mut_list[i].n++;
//				if (mut_list[i].n>mut_list[i].nmax){
	//				mut_list[i].list = (unsigned int *) realloc (mut_list[i].list, (mut_list[i].nmax+20)*sizeof (unsigned int)) ;
		//			mut_list[i].nmax += 20;
			//	}
				mut_list[i].list[mut_list[i].n-1] = total_mut;
				mut_list[i].list_hash[mut_list[i].n-1] = hash(total_mut);
				mut_list[i].two_copies[mut_list[i].n-1] = false;

                // (mut_list[i].list)[(mut_list[i].n)-1] = total_mut;
                // *(mut_list[i].list+mut_list[i].n-1) = total_mut;
			}

		}

	
		// print every 10th generation, starting at 9
		if (gen % 10 == 9){

			// every 1000th generation, starting at 9, open a new successes file and a new failures file 
            if (gen % 1000 == 9){

				// open file to print the current generation, the number of failures in this generation and the total number of mutations so far
	            sprintf(fffname, "sympatric_assortative_mating_NCH%dCONT_PROB%fNGEN%dMUT_PROB%ftheta%fseed_for_rand%dgen%d_number_of_failures.txt",NCH, CONT_PROB, NGEN, MUT_PROB, theta, seed_for_rand, gen); 
	            fff = fopen(fffname, "w");

	            // open file to print number of successes out of 1000 trials for couples of individuals in a list of distances
	            // (1, 10, 50, 100, 500)
                sprintf(ffname, "sympatric_assortative_mating_NCH%dCONT_PROB%fNGEN%dMUT_PROB%ftheta%fseed_for_rand%dgen%d_number_of_successes.txt",NCH, CONT_PROB, NGEN, MUT_PROB, theta, seed_for_rand, gen); 
	            ff = fopen(ffname, "w");

			}

			// print the current generation, the number of failures in this generation and the total number of mutations so far to a file
			fprintf (fff, "%d %d %d\n", gen, failures, total_mut);
						
			printf ("gen %d:\n", gen);   // print number of current generation

           	// print number of failures (number of inviable individuals produced) in current generation
			printf ("%d failures \n", failures);
			
			for (t=0; t<5; t++){
				success = 0;
				for (r=0; r<1000; r++){
                    
					// location of first parent
					unsigned int loc = ((unsigned int) rand()*(RAND_MAX) + (unsigned int) rand()) % NCH;
										
					// location of second parent
					int loc1 = (loc + try_dists[t]) % NCH;
										
					alive = true;
					success+= is_success (mut_list+loc, mut_list+loc1, &inv_cause1, &inv_cause2);
				}

				// print number of successes out of 1000 trials for couples of individuals in a list of distances
			    // (1, 10, 50, 100, 500) to a file
				//fprintf (ff, "%d ", success);
				Itoa(success, success_string, 10);
                fputs (success_string, ff);
				fputc (' ', ff);
				
			    printf ("%d: %d/1000\n", try_dists[t], success);
			}

			fputc ('\n', ff);
			
		    //printf ("times: %d %d %d %d tot mut: %d\n", t1_tot, t2_tot, t3_tot, t4_tot, total_mut);
			//printf ("tot1 %d tot11 %d tot2 %d tot3 %d tot4 %d\n",tot1,tot11,tot2,tot3,tot4);
		}
		//tot4 = tot4+clock()-t4;

        // close successes file and failures file 
        if ((gen % 1000 == 8) && (gen > 8)){
			fclose (ff);
	        fclose (fff);
		}

		
		// every 1000 generations print files containing information about successes/failures to reproduce and mutation
		
		if (gen % 1000 == 999){ 
			


			// print successes to reproduce as '1' and failures as '0' for couples of individuals in a list of distances to a file
			
			sprintf(fname, "sympatric_assortative_mating_NCH%dCONT_PROB%fNGEN%dMUT_PROB%ftheta%fseed_for_rand%dgen%d_just_0_1.txt",NCH, CONT_PROB, NGEN, MUT_PROB, theta, seed_for_rand, gen); 
			f = fopen(fname, "w");
			for (t=0; t<5; t++){
				for (r=0; r<NCH; r++){
					int loc1 = (r + try_dists[t]) % NCH;


					// the former command: fputc((is_success (mut_list+i, mut_list+loc1)==1)? '1':'0', f);

					if ((is_success (mut_list+r, mut_list+loc1, &inv_cause1, &inv_cause2) == 1)){
						fputc ('1', f);
						fputc (' ', f);
					}
					else{
						fputc ('0', f);
						fputc (' ', f);
					}
				}
				fprintf (f,"\n");
			}
			fclose (f);
			

			// for the same couples, for each failure, print the (first) couple of mutations whose combination caused inviability
			// or two 0's if there was a success (a viable offspring)
					
			sprintf(fname, "sympatric_assortative_mating_NCH%dCONT_PROB%fNGEN%dMUT_PROB%ftheta%fseed_for_rand%dgen%d_causative.txt",NCH, CONT_PROB, NGEN, MUT_PROB, theta, seed_for_rand, gen); 
			f = fopen(fname, "w");
			for (t=0; t<5; t++){
				for (r=0; r<NCH; r++){
					int loc1 = (r + try_dists[t]) % NCH;
					
					// the former command: fputc((is_success (mut_list+i, mut_list+loc1)==1)? '1':'0', f);

					if ((is_success (mut_list+r, mut_list+loc1, &inv_cause1, &inv_cause2) != 1)){
					
						// convert the mutation numbers into strings in base 10
					    Itoa (inv_cause1, mut_string1, 10);
						Itoa (inv_cause2, mut_string2, 10);

						// write the causative mutation numbers into the file
						fputs (mut_string1, f);
						fputc (' ', f);
						fputs (mut_string2, f);
						fputc (' ', f);
					}
					else{

						// write two 0's into the file
						fputs ("0 0 ", f);
					}

				}
				fprintf (f,"\n");
			}
			fclose (f);

			//sprintf(fname, "sympatric_assortative_mating_NCH%dCONT_PROB%fNGEN%dMUT_PROB%ftheta%fseed_for_rand%dgen%d_times.txt",NCH, CONT_PROB, NGEN, MUT_PROB, theta, seed_for_rand, gen); 
			//f = fopen(fname, "w");
			//fprintf (f,"tot1 %d tot2 %d tot3 %d tot4 %d\n",tot1,tot2,tot3,tot4);
			//fclose (f);




            // print the list of all the mutations of every individual
			
			// name of the file
			sprintf(fname_all_mutations, "sympatric_assortative_mating_NCH%dCONT_PROB%fNGEN%dMUT_PROB%ftheta%fseed_for_rand%dgen%d_all_mutations.txt",NCH, CONT_PROB, NGEN, MUT_PROB, theta, seed_for_rand, gen); 
			
			f = fopen(fname_all_mutations, "w");
			
			// go over all the individuals
			for (r=0; r<NCH; r++){

				// the number of mutations individual r has
				number_of_mutations = mut_list[r].n;

				// go over all of individual r's mutations
				for (int mutation_number=0; mutation_number<number_of_mutations; mutation_number++){

					// convert the mutation from unsigned int to string - not needed anymore because I use fprintf
					//convert(mut_list[i].list[mutation_number], this_mutation_string);
				
				    // Write the mutation into the file, and put a ' ' after it.
					fprintf (f, "%d ", mut_list[r].list[mutation_number]);
					//fputs (this_mutation_string, f);
					//fputc (' ', f);
					
				}

				// end the line before moving to the next individual
				fputc ('\n', f);

			}

			// close the file
			fclose (f);
			

			
            // for each of these mutations, print 0 if the individual has one copy of the mutation, and 1 if it has two copies.
			
			// name of the file
			sprintf(fname_all_mutations, "sympatric_assortative_mating_NCH%dCONT_PROB%fNGEN%dMUT_PROB%ftheta%fseed_for_rand%dgen%d_two_copies.txt",NCH, CONT_PROB, NGEN, MUT_PROB, theta, seed_for_rand, gen); 
			
			f = fopen(fname_all_mutations, "w");
			
			// go over all the individuals
			for (r=0; r<NCH; r++){

				// the number of mutations individual r has
				number_of_mutations = mut_list[r].n;

				// go over all of individual r's mutations
				for (int mutation_number=0; mutation_number<number_of_mutations; mutation_number++){

					if (mut_list[r].two_copies[mutation_number]){
						fputc ('1', f);
						fputc (' ', f);
					}
					else{
						fputc ('0', f);
						fputc (' ', f);
					}

				}

				// end the line before moving to the next individual
				fputc ('\n', f);

			}

			// close the file
			fclose (f);

		}

	}			
	
	
	
	return 0;
}