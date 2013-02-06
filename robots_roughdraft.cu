/*

robots.cu
Written by Spencer B Liberto
12 Dec 2012

An implementation of the simulator described in this paper:
http://stephane.magnenat.net/data/Evolutionary%20Conditions%20for%20the%20Emergence%20of%20Communication%20in%20Robots%20-%20Dario%20Floreano,%20Sara%20Mitri,%20St%C3%A9phane%20Magnenat,%20Laurent%20Keller%20-%20Current%20Biology%20-%202007.pdf

*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define NUMCLANS 10
#define ROBOTSPERCLAN 10
#define ARENASIDELENGTH 300
#define NUMGENERATIONS 100
#define VISION 100
#define FOODXY 100
#define POISINXY 200
#define FEEDINGDISTANCE 25
#define ROBOTRADIUS 5
#define NUMGENERATIONS 100
#define NUMCYCLESPERGEN 500
#define DATALOGFILENAME "datalog.txt"

typedef struct linkedListNode
{
	linkedListNode* prev;
	linkedListNode* next;
	int correspondingRobot;
} edibleQueueNode;

typedef struct anotherLinkedListNode
{
	anotherLinkedListNode* next;
	int correspondingRobot;
	int score;
} scoreListNode;

int main(){
	runSimulation();
	return 0;
}

/* ================================== */
/* == THE CYCLES OF THE SIMULATION == */
/* ================================== */

/* runs the entire simulation from start to finish */
void runSimulation(){
	int numRobots = NUMCLANS * ROBOTSPERCLAN;
	int statesLength = numRobots * 5;
	int phenotypesLength = numRobots * 240;

	int* prevPhenotypes = malloc(numPhenotypeBits * sizeof(int));
	int* prevStates = malloc(statesLength * sizeof(int));
	int* finalScores = malloc(numRobots * sizeof(int));

	int* devicePrevPhen = 0;
	int* devicePrevStat = 0;
	int* deviceFinalScores = 0;
	cudaMalloc((void**)devicePrevPhen, (phenotypesLength * sizeof(int)));
	cudaMalloc((void**)devicePrevStat, (statesLength * sizeof(int)));
	cudaMalloc((void**)deviceFinalScores, (numRobots * sizeof(int)));

	randomPhenotypes(prevPhenotypes);

	int i;
	for(i=0; i<NUMGENERATIONS; i++){
		randomStates(prevStates);
		cudaMemcpy(devicePrevPhen, prevPhenotypes, phenotypesLength, cudaMemcpyHostToDevice);
		cudaMemcpy(devicePrevStat, prevStates, statesLength, cudaMemcpyHostToDevice);
		generationCycle<<NUMCLANS,NUMROBOTS>>(devicePrevPhen, devicePrevStat, deviceFinalScores);
		cudaMemcpy(finalScores, deviceFinalScores, numRobots, cudaMemcpyDeviceToHost);
		produceNewGeneration(prevPhenotypes, finalScores, numRobots, i);
	}

	free(prevPhenotypes);
	free(prevStates);
	cudaFree(devicePrevPhen);
	cudaFree(devicePrevStat);
}

/* executes one generarion */
__global__ void generationCycle(int* globalPrevPhen, int* globalPrevStat, int* globalFinalScores){
	
	int tid = threadId.x;
	int bid = blockId.x;
	__shared__ int sharedPrevPhenotypes[blockPhenotypesLength];
	__shared__ int sharedPrevStates[blockPhenotypesLength];
	__shared__ edibleQueueNode edibleQueues[4]; /* [foodStart, foodEnd, PoisStart, PoisEnd] */
	__shared__ int finalScores[ROBOTSPERCLAN];

	if(tid == 0){
		int blockPhenotypesLength = ROBOTSPERCLAN * 240;
		int blockStatesLength = ROBOTSPERCLAN * 5;
		int* globalPrevBlockPhenotypes = (blockPhenotypesLength * blockIdx.x) + globalPrevPhen;
		int* globalPrevBlockStates = (blockStatesLength * blockIdx.x) + globalPrevStat;
		cudaMemcpy(sharedPrevPhenotypes, globalPrevBlockPhenotypes, blockPhenotypesLength, cudaMemcpyDeviceToDevice);
		cudaMemcpy(sharedPrevStates, globalPrevBlockStates, blockStatesLength, cudaMemcpyDeviceToDevice);
	}
	__syncthreads();

	int inputs[10];
	int outputs[3];
	int states[5];
	int phenotype[240];
	int translatedPhenotype[30];
	cudaMemcpy(states, sharedPrevStates, 5, cudaMemcpyDeviceToDevice);
	cudaMemcpy(phenotype, sharedPrevPhenotypes, 240, cudaMemcpyDeviceToDevice);
	translatePhenotype(phenotype, translatedPhenotype);

	int i;
	for(i=0; i<NUMCYCLESPERGEN; i++)
		ioCycle(tid, states, sharedPrevStates, inputs, outputs, translatedPhenotype, edibleQueues);
	for(i=0; i<ROBOTSPERCLAN; i++)
		finalScores[i] = sharedPrevStates[(5*i)+4];
	if(tid=0)
		cudaMemcpy(globalFinalScores[(bid*ROBOTSPERCLAN)], finalScores, ROBOTSPERCLAN, cudaMemcpyDeviceToDevice);
	__syncthreads();
}

/* executes one iocycle of a generation */
__device__ void ioCycle(int tid, int states[5], int* sharedPrevStates, int inputs[10], int outputs[3], int translatedPhenotype[30], edibleQueueNode* edibleQueues){
	/* One IOCycle */
	int bstates[5];
	int i;
	for(i=0; i<ROBOTSPERCLAN; i++){
		if (i != tid){
			cudaMemcpy(bstates, (sharedPrevStates+(5*i)), 5, cudaMemcpyDeviceToDevice);
			updateInputsRobots(states, bstates, inputs);
		}
	}
	updateInputsEdibles(states, inputs, edibleQueues, tid);
	updateOutputs(outputs, inputs, translatedPhenotype, states);
	updateStatesXYOB(states, outputs);
	if(tid == 0){
		tallyPoints(sharedPrevStates, edibleQueues)	
	}
	__syncthreads();
}

/* ========================================== */
/* == PRE-GENERATION RANDOM INITIAL VALUES == */
/* ========================================== */

/* produces a list of random initial states */
void randomStates(int* initialStates){
	int numRobots = NUMCLANS * ROBOTSPERCLAN;
	int i;
	for(i=0; i<numRobots; i++){
		&(initialStates++) = random(ARENASIDELENGTH);
		&(initialStates++) = random(ARENASIDELENGTH);
		&(initialStates++) = random(8) * pi / 4;
		&(initialStates++) = random(11) / 10;
		&(initialStates++) = 0;
	}
}

/* produces a list of random binart phenotypes */
void randomPhenotypes(int* initialPhenotypes){
	int numPhenotypeBits = NUMCLANS * ROBOTSPERCLAN * 480;
	int i;
	for(i=0; i<totalPhenotypeBits; i++)
		&(initialPhenotypes++) = random(2);
}

/* ======================================== */
/* == POST-GENERATION MATING AND LOGGING == */
/* ======================================== */

/* finds top 20% scorers, mates them, and produces a new generation of phenotypes. Also prints to log. */
void produceNewGeneration(int* oldPhenotypes, int* finalScores, int numRobots, int* states, int generationNumber){
	int numMedalists = floor(numRobots/5);

	/* rank scores */
	scoreListNode* scoreList = (scoreListNode *) malloc(sizeof(scoreListNode));
	scoreList->next = NULL;
	scoreList->correspondingRobot = 0;
	scoreList->score = finalScores[0];
	scoreListNode* iterScoreList;
	/* For each score, traverse list from top until lower score found */
	int i, j;
	for(i=1; i<numRobots; i++){
		scoreListNode* newNode = (scoreListNode *) malloc(sizeof(scoreListNode));
		newNode->correspondingRobot = i;
		newNode->score = finalScores[i];
		iterScoreList = scoreList;
		j=0
		while((j<i) && (finalScores[i] < iterScoreList->score)){
			iterScoreList = iterScoreList->next;
		}
		newNode->next = iterScoreList->next;
		iterScoreList->next = newNode;
	}

	/*  Output data */
	printScoreList(scoreList, numRobots, oldPhenotypes, states, generationNumber);

	/* find top scorers ("medalists") */
	int medalists[numMedalists];
	iterScoreList = scoreList;
	for(i=0; i<numMedalists; i++){
		medalist[i] = iterScoreList->correspondingRobot;
		iterScoreList = iterScoreList->next;
	}

	/* mate top scorers */
	int* deviceOldPhen;
	int* deviceNewPhen;
	int* deviceMedalists;
	cudaMalloc((void**)deviceOldPhen, (numRobots * 240 * sizeof(int)));
	cudaMalloc((void**)deviceNewPhen, (numRobots * 240 * sizeof(int)));
	cudaMalloc((void**)deviceMedalists, (numRobots * sizeof(int)));
	cudaMemcpy(deviceOldPhen, oldPhenotypes, (240*numRobots), cudaMemcpyHostToDevice);
	cudaMemcpy(deviceMedalists, medalists, numMedalists, cudaMemcpyHostToDevice);
	mateRobotsRandomly<<1,numRobots>>(deviceOldPhen, deviceNewPhen, medalists, numMedalists);
	cudaMemcpy(oldPhenotypes, deviceNewPhen, (240*numRobots), cudaMemcpyDeviceToHost);
}

/* randomly mates specified robots */
__global__ void mateRobotsRandomly(int* oldPhen, int* newPhen, int* medalists, int numMedalists){

	/* Get information about parents */
	int tid = threadIdx.x;
	int* child = malloc(240*sizeof(int));
	int parenta = 0;
	int parentb = 0;
	while(parenta==parentb){
		parenta = medalist + (random(numMedalists));
		parentb = medalist + (random(numMedalists));
	}
	int* iterA = (oldPhen + 240*parenta);
	int* iterB = (oldPhen + 240*parentb);
	int* iterChild = child;
	cudaMalloc((void**)iterA, (numRobots * 240 * sizeof(int)));
	cudaMalloc((void**)iterB, (numRobots * 240 * sizeof(int)));
	cudaMemcpy(iterA, (oldPhen + (240*parenta)), 240, cudaMemcpyDeviceToDevice);
	cudaMemcpy(iterB, (oldPhen + (240*parentb)), 240, cudaMemcpyDeviceToDevice);

	/* Populate child */
	int i, currentA;
	for(i=0; i<240; i++){
		currentA = &iterA;
		if((&(iterA++)==(&(iterB++))) || (tid%2 == 0))
			&(child++) = currentA;
		else
			&(child++) = 1 - currentA;
	}

	/* Send child to global */
	cudaMemcpy((newPhen+(240*tid)), child, 240, cudaMemcpyDeviceToDevice);
}

/*  Outputs data about the current generation */
void printScoreList(scoreListNode* scoreList, int numRobots, int* phenotypes, int* states, int generationNumber){

	FILE *fp;
	fp=fopen(DATALOGFILENAME, "a");

	int translatedPhen[30];
	int currentRobot;
	int* iterPhen;
	int* iterStates;
	int i,j,k;
	fprintf(fp, "generation %d\n\nDNA\nx-coordinate y-coordinate orientation brightness final-score\n\n", generationNumber);
	for(i=0; i<numRobots; i++){
		currentRobot = scoreList->correspondingRobot;
		iterPhen = phenotypes + (240*currentRobot);
		iterStates = states + (5*currentRobot);
		for(j=0; j<30; j++){
			for(k=0; k<8; k++){
				fprintf(fp,"%d", &iterPhen);
				iterPhen++
			}
			fprintf(fp, " ");
		}
		printf("\n");
		for(j=0; j<4; j++){
			fprintf(fp,"%d   ", &(iterStates++));
		}
		pfrintf(fp,"%d\n", (scoreList->score);
		scoreList = scoreList->next;
	}

	fclose(fp);
}

/* ================================= */
/* == END-OF-CYCLE POINT COUNTING == */
/* ================================= */

/* adds or removes a robot from an edible-queue */
__device__ void queueChanging(int task, int robotInQuestion, edibleQueueNode* edibleQueues){

	switch(task){
		/* add to foodQueue */
		case 0:
			edibleQueueNode* newNode;
			cudaMalloc((void**)newNode, sizeof(edibleQueueNode));
			newNode->prev = edibleQueues[1];
			newNode->next = NULL;
			newNode->correspondingRobot = robotInQuestion;
			(edibleQueues[1])->next = newNode;
			edibleQueues[1] = newNode;
			break;
		/* remove from foodQueue */
		case 1:
			edibleQueueNode* currentPointer = edibleQueues[0];
			if(currentPointer != NULL){
				while((currentPointer->correspondingRobot != robotInQuestion) || (currentPointer->next != NULL)){
					currentPointer = currentPointer->next;
				}
				if(currentPointer->correspondingRobot == robotInQuestion){
					if(currentPointer->next != NULL)
							(currentPointer->next)->prev = currentPointer->prev;
					(currentPointer->prev)->next = currentPointer->next;
				}
			}
			break;
		/* add to foodQueue */
		case 2:
			edibleQueueNode* newNode;
			cudaMalloc((void**)newNode, sizeof(edibleQueueNode));
			newNode->prev = edibleQueues[3];
			newNode->next = NULL;
			newNode->correspondingRobot = robotInQuestion;
			(edibleQueues[3])->next = newNode;
			edibleQueues[3] = newNode;
			break;
		/* remove from foodQueue */
		case 3:
			edibleQueueNode* currentPointer = edibleQueues[1];
			if(currentPointer != NULL){
				while((currentPointer->correspondingRobot != robotInQuestion) || (currentPointer->next != NULL)){
					currentPointer = currentPointer->next;
				}
				if(currentPointer->correspondingRobot == robotInQuestion){
					if(currentPointer->next != NULL)
							(currentPointer->next)->prev = currentPointer->prev;
					(currentPointer->prev)->next = currentPointer->next;
				}
			}
			break;
		default:
			fprintf(stderr, "Error: queueChanging recieved an unknown task.");
	}
}

/* tallies the points after an iocycle */
__device__ void tallyPoints(int* states, edibleQueueNode* edibleQueues){
	edibleQueueNode* currentPointer = edibleQueues[0];
	
	int i;
	for(i=0; i<5; i++){
		if(currentPointer!=NULL)
			states[(((currentPointer->correspondingRobot)*5)+4)] += 1;
		else
			break;
	}

	edibleQueueNode* currentPointer = edibleQueues[2];
	for(i=0; i<5; i++){
		if(currentPointer!=NULL)
			states[(((currentPointer->correspondingRobot)*5)+4)] -= 1;
		else
			break;
	}
}

/* ========================== */
/* == IOCYCLE CALCULATIONS == */
/* ========================== */

/* updates outputs for one robot, based upon inputs and phenotype */
__device__ void updateOutputs(int outputs[3], int inputs[10], int phenotype[30], int states[5]){
	for(iterOutputs=0; iterOutputs<3; iterOutputs++){
		for(iterInputs=0; iterInputs<10; iterInputs++)
			outputs[iterOutputs] += tanh(inputs[iterInputs] * phenotype[(10 * iterOutputs) + iterInputs]);
		outputs[iterOutputs] = tanh(outputs[iterOutputs]);
	}
	states[3] = outputs;
}

/* updates inputs for one robot, based upon one other robots outputs */
__device__ void updateInputsRobots(int states[5], int bstates[5], int inputs[10]){
	int xdiff, ydiff, distance, inputBrightness, theta, bucket;

	/* differences in x and y axis */
	xdiff = states[0] - bstates[0];
	ydiff = states[1] - bstates[1];

	/* Calculate percieved brightness, break if <0 */
	distance = sqrt(pow(xdiff,2) + pow(ydiff,2));
	inputBrightness = states[3] - (distance/VISION);
	if inputBrightness<0 {
		break;
		printf("Break Error: updateInputsRelative\n")
	}

	/* calculate which bucket the light will be picked up by */
	if xdiff==0{
		if ydiff>0
			theta = pi/2;
		else 
			theta = 3*pi/2;	
	} else {
		theta = atan(ydiff/xdiff);
		if xdiff<0
			theta += pi;
	}
	bucket = floor(((states[2]+theta) *4 /pi) %8);

	/* add light value to appropriate bucket */
	inputs[bucket] += inputBrightness;
}

/* updates inputs, based upon edibles outputs */
__device__ void updateInputsEdibles(int states[5], int inputs[10], edibleQueueNode* edibleQueues, int robotInQuestion){
	int xdiff, ydiff, distance, inputBrightness, theta, bucket;

	/* CALCULATE INPUTS FROM FOOD */

	xdiff = states[0] - FOODXY;
	ydiff = states[1] - FOODXY;
	distance = sqrt(pow(xdiff,2) + pow(ydiff,2));

	/* check if on food */
	if((distance<=FEEDINGDISTANCE) && (states[8]!=1))
		queueChanging(0, robotInQuestion, edibleQueues);
	else {
		if(states[8]==1)
			queueChanging(1, robotInQuestion, edibleQueues);
		states[8] *= 0.95;
	}

	/* calculate percieved brightness from food */
	inputBrightness = 1 - (distance/(2*ARENASIDELENGTH));

	/* calculate bucket which will recieve the light input */
	if xdiff==0{
		if ydiff>0
			theta = pi/2;
		else
			theta = 3*pi/2;
	} else {
		theta = atan(ydiff/xdiff);
		if xdiff<0
			theta += pi; 
	}
	bucket = ((states[2]+theta) *4 /pi) %8;

	/* Subtract input value from appropriate light bucket */
	inputs[bucket] -= inputBrightness;

	/* CALCULATE INPUTS FROM POISIN */

	xdiff = states[0] - POISINXY;
	ydiff = states[1] - POISINXY;
	distance = sqrt(pow(xdiff,2) + pow(ydiff,2));

	/* check if on poisin */
	if((distance<=FEEDINGDISTANCE) && (states[9]!=1))
		queueChanging(2, robotInQuestion, edibleQueues);
	else {
		if(states[9]==1)
			queueChanging(3, robotInQuestion, edibleQueues);
		states[9] *= 0.95;
	}

	/* calculate percieved brightness from poisin */
	inputBrightness = 1 - (distance/(2*ARENASIDELENGTH));

	if xdiff==0{
		if ydiff>0
			theta = pi/2;
		else
			theta = 3*pi/2;
	} else {
		theta = atan(ydiff/xdiff);
		if xdiff<0
			theta+=pi;
	}
	bucket = floor(((states[2]+theta) *4 /pi) %8);

	/* Subtract input value from appropriate light bucket */
	inputs[bucket] -= inputBrightness;
}

/* translates a binary phenotype into decimal */
__device__ void translatePhenotype(int untranslated[240], int translated[30]){
	int* untransIter = untranslated[0];
	int* transIter = translated[0];

	int i, j;
	for(i=0; i<30; i++){
		*transIter = 0;
		for(j=0; j<8; j++){
			if(*(untransIter++))
				*transIter += pow(2, (7-j));
		}
		transIter++;
	}
}

/* updates the states of a robot, based upon its outputs */
__device__ void updateStatesXYOB(int states[5], int outputs[3]){
	int arcRadius, dir, relativeAngle, m, n, xmove, ymove, newx, newy, newo;

	/* Calculate new x and y coordinates
		This is a quick fix, I would love an algorithm written by
		someone who actually knows trigonometry */
	if (outputs[0] == outputs[1]) {
		if outputs[0] == 0{
			break;
			printf("Break Error: updateStatesXYO");
		}
		arcRadius = pi * ROBOTRADIUS;
		relativeAngle = 0;
	} else {
		if (outputs[0] > outputs[1]) {
			m = outputs[0];
			n = outputs[1];
			dir = 0;
		} else {
			m = outputs[1];
			n = outputs[0];
			dir = 1;
		}
		relativeAngle = (m-n)*pi;
		if n==0
			arcRadius = ROBOTRADIUS;
		else
			arcRadius = ROBOTRADIUS*m/n;
	}	
	absoluteAngle = (relativeAngle + states[2]) % (2*pi);
	if (dir)
		absoluteAngle *= -1;
	xmove = sin(absoluteAngle) * arcRadius;
	ymove = (cos(absoluteAngle) - 1) * arcRadius;
	newx = xmove + states[0];
	newy = ymove + states[1];

	/* Robots shouldn't be allowed out of the boundaries.
		This boundry-collision algorithm is not very physics-like,
		I would love an algorithm that could be more accurate */
	if (newx <= 0)
		newx = 0;
	if (newx >= ARENASIDELENGTH)
		newx = ARENASIDELENGTH;
	if (newy <= 0)
		newy = 0;
	if (newy >= ARENASIDELENGTH)
		newy = ARENASIDELENGTH;

	/* Calculate new orientation */
	if (dir)

		calc = (pi*2);

		newo = (absoluteAngle - (pi/2)) % calc;
	else
		newo = (absoluteAngle + (pi/2)) % calc;

	/* Assign new values */
	states[0] = newx;
	states[1] = newy;
	states[2] = newo;
	states[3] = outputs[2]; /* The brightness can directly transfer */
}