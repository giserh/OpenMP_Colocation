#pragma warning(disable : 4996)
#include "includes.h"
#include<omp.h>
#include "colocationFinder.h"
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/compare.hpp>
#include <string.h>
#include <algorithm>  
#include <iterator>  
#include<valarray>
#include<vector>
#include<numeric>
#include<functional>
#include<unordered_map>


using namespace boost;
using namespace std;

colocationFinder::colocationFinder(void) {
	totalCandidatePatterns = 0;
	totalFilteredPatterns = 0;
	totalFilterTime = 0.0;
	totalRefineTime = 0.0;
	degree2FilterTime = 0.0;
	totalKernelLaunches = 0;
	h_coloc = new Integer[10];
	instanceMax = 30000000;
	table1Max = instanceMax * 5;
	intermediateMax = 120000000 * 10;
	h_Indexes = (Index64*)malloc(instanceMax*sizeof(Index64));
	preprocessingTime = 0.0;
	ioTime = 0.0;
}
colocationFinder::~colocationFinder(void) {
	featureTypes.clear();
	featureTypes.shrink_to_fit();
	free(instanceList);
	delete[] featureInstanceStart;
	delete[] featureInstanceEnd;
	delete[] featureInstanceCount;
	free(instanceLocationX);
	free(instanceLocationY);
	h_candiColocations.clear();
	h_candiColocations.shrink_to_fit();
	
}
void colocationFinder::LoadParameter(std::string parameterFileName) {
	std::string FileName = /*location+ */parameterFileName;
	std::ifstream ifs(FileName.c_str());
	if (false == ifs.is_open())
	{
		std::cout << "Parameter file missing..." << std::endl;
		exit(0);
		//return CCT_FILEERR;
	}
	ifs >> parameter.thresholdDistance;
	ifs >> parameter.PIthreshold;
	ifs >> parameter.FilterON_OFF;
	ifs.close();
	parameter.squaredDistanceThreshold = parameter.thresholdDistance* parameter.thresholdDistance;
}
void colocationFinder::degree2CandidateColocationGenerator() {
	for (Integer i = 0; i < featureTypes.size() - 1; i++) {
		for (Integer j = i + 1; j < featureTypes.size(); j++) {
			h_candiColocations.push_back(i);
			h_candiColocations.push_back(j);
		}
	}
}
Integer mode(Integer a[], Integer n) {
	Integer maxValue = 0, maxCount = 0, i, j;

	for (i = 0; i < n; ++i) {
		Integer count = 0;

		for (j = 0; j < n; ++j) {
			if (a[j] == a[i])
				++count;
		}

		if (count > maxCount) {
			maxCount = count;
			maxValue = a[i];
		}
	}

	return maxValue;
}
void colocationFinder::setGrid(Real maxX, Real maxY, Real minX, Real minY) {
	gridStructure.computeZone.minBound.x = minX;
	gridStructure.computeZone.minBound.y = minY;
	gridStructure.computeZone.maxBound.x = maxX;
	gridStructure.computeZone.maxBound.y = maxY + parameter.thresholdDistance;
	gridStructure.cellWidth = parameter.thresholdDistance;
	scalar2 gs;
	gs.x = (gridStructure.computeZone.maxBound.x - gridStructure.computeZone.minBound.x) / gridStructure.cellWidth;
	gs.y = (gridStructure.computeZone.maxBound.y - gridStructure.computeZone.minBound.y) / gridStructure.cellWidth;
	gridStructure.gridSize.x = (Integer)ceil(gs.x);
	gridStructure.gridSize.y = (Integer)ceil(gs.y);
	gridStructure.totalCells = gridStructure.gridSize.x * gridStructure.gridSize.y;
	Integer cSize = gridStructure.totalCells*maxFeaturesNum;
	std::cout << "cSize" << cSize;
 	cellEventInstanceCount = (SInteger*)malloc(cSize*sizeof(SInteger));
	instanceCellIDs = (Integer*)malloc(maxInstancesNum*sizeof(Integer));
	cellBasedSortedIndex = (Integer*)malloc(maxInstancesNum*sizeof(Integer));
	
	for (size_t i = 0; i < gridStructure.totalCells*maxFeaturesNum; i++) {
		cellEventInstanceCount[i] = 0;
	}
	vector< pair<Integer, Integer> >v;
	for (size_t i = 0; i < maxInstancesNum; i++) {
		Integer typeID = instanceList[i];
		Integer cellID = getCellID(instanceLocationX[i], instanceLocationY[i]);
		instanceCellIDs[i] = cellID;
		Integer value = cellEventInstanceCount[cellID*maxFeaturesNum + typeID];
		value++;
		cellEventInstanceCount[cellID*maxFeaturesNum + typeID] = value;
		cellBasedSortedIndex[i] = i;
		v.push_back(make_pair(cellID, i));
	}

	std::sort(std::begin(v), std::end(v));
	Integer cnt = 0;
	for (Integer t = 0; t < v.size(); t++) {
		cellBasedSortedIndex[cnt] = v[t].second;
		cnt++;
	}
	v.clear();
	v.shrink_to_fit();
	size_t posCounter = 0;
	size_t maxEnd = 0;
	cellEventInstanceStart = (SInteger*)malloc(cSize*sizeof(SInteger));
	cellEventInstanceEnd = (SInteger*)malloc(cSize*sizeof(SInteger));
	for (size_t CID = 0; CID < gridStructure.totalCells; CID++) {
		for (size_t FID = 0; FID < maxFeaturesNum; FID++) {
			Integer countValue = cellEventInstanceCount[CID*maxFeaturesNum + FID];
			if (countValue > 0) {
				cellEventInstanceStart[CID*maxFeaturesNum + FID] = posCounter;
				cellEventInstanceEnd[CID*maxFeaturesNum + FID] = posCounter + countValue - 1;
				maxEnd = posCounter + countValue - 1;
				posCounter = posCounter + countValue;
			}
			else {
				cellEventInstanceStart[CID*maxFeaturesNum + FID] = -1;
				cellEventInstanceEnd[CID*maxFeaturesNum + FID] = -1;
			}
		}
	}
}
Integer colocationFinder::getCellID(Real x, Real y) {
	Real RelativeX = x - gridStructure.computeZone.minBound.x;
	Real RelativeY = y - gridStructure.computeZone.minBound.y;

	RelativeX /= gridStructure.cellWidth;
	RelativeY /= gridStructure.cellWidth;
	Integer i = (Integer)RelativeX;
	Integer j = (Integer)RelativeY;
	const Integer CellID = (gridStructure.gridSize.x * j) + i;
	return CellID;
}

void colocationFinder::populateData(std::string datasetFilename) { //pass data.csv as argument

	auto begin = get_time::now();															   //variable decleration and defination
	vector<string> vec;
	string line;
	size_t entryCounter = 0;
	size_t eventTypeCounter = 0;
	size_t instanceCounter = 0;
	size_t lastEventID = -1;
	vector<struct sFeatureStats> featureStats;
	Real maxX, maxY, minX, minY;
	string data(/*location +*/ datasetFilename);

	//File Reading
	ifstream in(data.c_str());
	typedef tokenizer< escaped_list_separator<char> > Tokenizer;
	getline(in, line); // for instance count
	maxInstancesNum = stoi(line);
	instanceLocationX = (Real*)malloc(maxInstancesNum*sizeof(Real));
	instanceLocationY = (Real*)malloc(maxInstancesNum*sizeof(Real));
	instanceList = (Integer*)malloc(maxInstancesNum*sizeof(Integer));
	getline(in, line);
	while (!line.empty())
	{
		Tokenizer tok(line);
		vec.assign(tok.begin(), tok.end());
		if (lastEventID == -1) {
			lastEventID = stoi(vec[0]);
			featureTypes.push_back(vec[0]);
			eventTypeCounter++;
		}
		else if (stoi(vec[0]) != lastEventID) {
			struct sFeatureStats tempStat;
			tempStat.start = entryCounter - instanceCounter;
			tempStat.end = entryCounter - 1;
			tempStat.count = instanceCounter;
			featureStats.push_back(tempStat);
			lastEventID = stoi(vec[0]);
			instanceCounter = 0;
			featureTypes.push_back(vec[0]);
			eventTypeCounter++;
		}
		instanceList[entryCounter] = eventTypeCounter - 1;  //hold the feature id of each instances... all the instances are sorted by type before read.
		instanceLocationX[entryCounter] = stof(vec[1]);
		instanceLocationY[entryCounter] = stof(vec[2]);
		entryCounter++;
		instanceCounter++;
		getline(in, line);
	}
	getline(in, line);
	Tokenizer tok(line);
	vec.assign(tok.begin(), tok.end());
	minX = stof(vec[0]);
	minY = stof(vec[1]);
	maxX = stof(vec[2]);
	maxY = stof(vec[3]);
	//for last feature
	struct sFeatureStats tempStat;
	tempStat.start = entryCounter - instanceCounter;
	tempStat.end = entryCounter - 1;
	tempStat.count = instanceCounter;
	featureStats.push_back(tempStat);

	//converting vector to array
	Integer featureSize = featureStats.size();
	featureInstanceStart = new Integer[featureSize];
	featureInstanceEnd = new Integer[featureSize];
	featureInstanceCount = new Integer[featureSize];
	Integer MaxCount= 0; 
	for (size_t i = 0; i < featureStats.size(); i++) {
		featureInstanceStart[i] = featureStats[i].start;
		featureInstanceEnd[i] = featureStats[i].end;
		featureInstanceCount[i] = featureStats[i].count;
		if (featureInstanceCount[i] > MaxCount) {
			MaxCount = featureInstanceCount[i];
		}
	}
	mapMax = MaxCount * 8;
	maxFeaturesNum = featureSize;
	in.close();
	auto last = get_time::now();
	auto difference = last - begin;
	Real msecs = chrono::duration_cast<chrono::milliseconds>(difference).count();
	ioTime += msecs;
	begin = get_time::now();
	std::cout << std::endl << "Data Load time (PopulateData Function) : " << msecs << " milliseconds " << endl;
	setGrid(maxX, maxY, minX, minY);
	last = get_time::now();
	difference = last - begin;
	msecs = chrono::duration_cast<chrono::milliseconds>(difference).count();
	preprocessingTime += msecs;
	std::cout << std::endl << "Grid Processing time : " << msecs << " milliseconds " << endl;
}

Real colocationFinder::DistanceInMeters(size_t idx1, size_t idx2) {
	Real x2 = instanceLocationX[idx2];
	Real x1 = instanceLocationX[idx1];
	Real y2 = instanceLocationY[idx2];
	Real y1 = instanceLocationY[idx1];
	Real squareDistance = pow((x2 - x1), 2) + pow((y2 - y1), 2);
	return squareDistance;
}

SInteger colocationFinder::getNeighborCellID(Integer x, Integer y, Integer TID) {
	SInteger NCID;
	NCID = y*gridStructure.gridSize.x + x;
	if (TID%gridStructure.gridSize.x == 0) {
		if ((NCID + 1) % gridStructure.gridSize.x == 0) {
			return -1;
		}
		return NCID;
	}
	else if ((TID + 1) % gridStructure.gridSize.x == 0) {
		if (NCID % gridStructure.gridSize.x == 0) {
			return -1;
		}
		return NCID;
	}
	else {
		return NCID;
	}
}
void colocationFinder::degree2Processing() {
	//std::cout << "\n\n\nFor degree 2......" << std::endl;
	auto start = get_time::now();
	//timeLog(2, "candidate colocgen", "start", 0);
	candiColocCounter = 0;
	for (size_t i = 1; i <= featureTypes.size() - 1; ++i)
	{
		candiColocCounter += i;
	}
	Integer degree = 2;
	totalCandidatePatterns += candiColocCounter;
	degree2CandidateColocationGenerator();
	auto end = get_time::now();
	auto diff = end - start;
	h_prevalentSlotCounts = (Integer**)malloc((candiColocCounter)*sizeof(Integer*));
	estimatePrevalent = candiColocCounter;
	slotLimit = instanceMax;
	for (int i = 0; i < estimatePrevalent; i++) {
		h_prevalentSlotCounts[i] = (Integer*)malloc(slotLimit*sizeof(Integer));
	}
	//timeLog(2, "candidate colocgen " + std::to_string(candiColocCounter), "end", chrono::duration_cast<sec>(diff).count());
	//std::cout << "degree2CandidateColocationGenerator() took  " << chrono::duration_cast<sec>(diff).count() << " seconds " << endl;
	//***********CUDA Implementation: instance table generation*************
	StatusType status;

	h_prevelentColocationCount = 0;
	//timeLog(2, "multiResolution pruning + Instance Table Generation ", "start", 0);
	Integer filterpruneCounter = 0;
	//compLog(degree, "generate Table Instance Start with total candidate Colocation patterns = " + std::to_string(candiColocCounter));
	for (size_t i = 0; i < candiColocCounter; i++) {
		//std::cout << "Degree " << degree << ": Candidate pattern " << i + 1 << "/" << candiColocCounter << " start" << std::endl;
		size_t firstFeatureIndex = h_candiColocations[i*degree];
		size_t secondFeatureIndex = h_candiColocations[i*degree + 1];
		size_t start = featureInstanceStart[firstFeatureIndex];
		size_t end = featureInstanceEnd[firstFeatureIndex];
		size_t secondstart = featureInstanceStart[secondFeatureIndex];
		for (Integer t = 0; t < degree; t++) {
			h_coloc[t] = h_candiColocations[i*degree + t];
		}
		//Filtering Mechanism
		if (parameter.FilterON_OFF == 1) {
			auto fstart = get_time::now();
			Real upperbound = filterCPP(h_coloc, degree);
			if (upperbound < parameter.PIthreshold) {
				filterpruneCounter++;
				//compLog(i + 1, degree, candiColocCounter, "pruned", " upperbound: ", upperbound);
				//std::cout << "Degree " << degree << ": Candidate pattern " << i + 1 << "/" << candiColocCounter << " end (filter based pruned)" << std::endl;
				auto fend = get_time::now();
				auto fdiff = fend - fstart;
				Real msecs = chrono::duration_cast<chrono::milliseconds>(fdiff).count();
				totalFilterTime += (msecs / 1000.0);
				continue;
			}
			//compLog(i, degree, candiColocCounter, "", " upperbound: ", upperbound);
			auto fend = get_time::now();
			auto fdiff = fend - fstart;
			Real msecs = chrono::duration_cast<chrono::milliseconds>(fdiff).count();
			totalFilterTime += (msecs / 1000.0);
		}
		//Filtering Mechanism Ends
		auto rStart = get_time::now();
		size_t totalInstances = 0;
		size_t table1InstanceCount = end - start + 1;
		size_t bitmapSize = featureInstanceCount[firstFeatureIndex] + featureInstanceCount[secondFeatureIndex];
		vector<Integer> bitmap(bitmapSize, 0);
		if (table1InstanceCount > slotLimit) {
			free(h_prevalentSlotCounts[h_prevelentColocationCount]);
			h_prevalentSlotCounts[h_prevelentColocationCount] = (Integer*)malloc(table1InstanceCount*sizeof(Integer));
			slotLimit = table1InstanceCount;
		}
		memset(h_prevalentSlotCounts[h_prevelentColocationCount], 0, sizeof(Integer)*table1InstanceCount);
		degree2TableSlotCounter_Kernel(table1InstanceCount, h_prevalentSlotCounts[h_prevelentColocationCount], start, end, secondstart, bitmap, secondFeatureIndex);
		totalKernelLaunches++;
		Real PI = getPI(bitmap, degree, h_coloc);
		if (PI < parameter.PIthreshold) {
				continue;
		}
		h_prevalantColocations.push_back(firstFeatureIndex);
		h_prevalantColocations.push_back(secondFeatureIndex);
		h_prevelentColocationCount++;
		auto rend = get_time::now();
		auto rdiff = rend - rStart;
		Real mrsecs = chrono::duration_cast<chrono::milliseconds>(rdiff).count();
		totalRefineTime += (mrsecs / 1000.0);
	}
	hasInstanceTable.clear();
	hasInstanceTable.shrink_to_fit();
	for (Integer i = 0; i < h_prevelentColocationCount; i++) {
		hasInstanceTable.push_back(0);
	}
	candidateColocationGeneral(degree+1);
	if (candiColocCounter == 0) {
		needInstanceTable = false;
		degree2FilterTime = totalFilterTime;
		totalFilteredPatterns += filterpruneCounter;
		return;
	}
	Integer degkplus1 = degree + 1;
	for (Integer i = 0; i < candiColocCounter; i++) {
		vector<Integer> combiColoc1;
		for (size_t t = 0; t < degkplus1; t++) {
			if (t < degkplus1 - 2) {
				combiColoc1.push_back(h_candiColocations[i*degkplus1 + t]);
			}
		}
		combiColoc1.push_back(h_candiColocations[i*degkplus1 + (degkplus1 - 2)]);

		Integer table1Index = getIndex(combiColoc1, degkplus1 - 1);
		hasInstanceTable[table1Index] = 1;
	}
	h_prevalentInstanceTable = (Integer**)malloc(h_prevelentColocationCount*sizeof(Integer*));
	auto refineStart = get_time::now();
	for (size_t i = 0; i < h_prevelentColocationCount; i++) {
		if (hasInstanceTable[i] == 1) {
			size_t firstFeatureIndex = h_prevalantColocations[i*degree];
			size_t secondFeatureIndex = h_prevalantColocations[i*degree + 1];
			size_t start = featureInstanceStart[firstFeatureIndex];
			size_t end = featureInstanceEnd[firstFeatureIndex];
			size_t secondstart = featureInstanceStart[secondFeatureIndex];
			size_t table1InstanceCount = end - start + 1;
			if (table1InstanceCount > instanceMax) {
				free(h_Indexes);
				h_Indexes = (Index64*)malloc(instanceMax*sizeof(Index64));
			}
			exclusive_scan(i, table1InstanceCount);
			Integer totalInstances = getTotalInstances(h_prevalentSlotCounts[i], table1InstanceCount);
			//std::cout << "Total Instance size =" << totalInstances << std::endl;
			size_t totalSize = totalInstances*degree;
			#pragma omp parallel for
			for (int TID = 0; TID < table1InstanceCount; TID++) {
				size_t temp = h_Indexes[TID];
				size_t value = temp * degree;
				h_Indexes[TID] = value;
			}
			h_prevalentInstanceTable[i] = (Integer*)malloc(totalInstances*degree*sizeof(Integer));
			degree2TableGenerator_Kernel(table1InstanceCount,h_prevalentSlotCounts[i], h_Indexes, h_prevalentInstanceTable[i], start, end, secondFeatureIndex);
			totalKernelLaunches++;

			h_prevalentInstanceTableSize.push_back(totalInstances);
		}
		else {
			h_prevalentInstanceTableSize.push_back(0);
		}
	}
	auto refineEnd = get_time::now();
	auto refineDiff = refineEnd - refineStart;
	Real mrsecs = chrono::duration_cast<chrono::milliseconds>(refineDiff).count();
	totalRefineTime += (mrsecs / 1000.0);

	degree2FilterTime = totalFilterTime;
	totalFilteredPatterns += filterpruneCounter;
}


void colocationFinder::tableGenRequired(Integer degree) {
	Integer lastIndex = degree - 2;
	bool flag;
	needInstanceTable = false;
	for (size_t i = 0; i < candiColocCounter; i++) {
		flag = true;
		for (size_t j = i + 1; j < candiColocCounter; j++) {
			for (size_t k = 0; k < lastIndex; k++) {
				if (h_candiColocations[i*(degree - 1) + k] != h_candiColocations[j*(degree - 1) + k]) {
					flag = false;
					break;
				}
			}
			if (flag && h_candiColocations[i*(degree - 1) + lastIndex]< h_candiColocations[j*(degree - 1) + lastIndex]) {
				vector<Integer> inter;
				for (Integer t = 0; t < degree - 1; t++) {
					inter.push_back(h_candiColocations[i*(degree - 1) + t]);
				}
				inter.push_back(h_candiColocations[j*(degree - 1) + lastIndex]);
				//module to generate k-1 degree subset
				flag = generateandCheckSubsets2(inter, degree);
				clearSubsetVectors();
				if (!flag) {
					break;
				}
				else {
					needInstanceTable = true;
					return;
				}
			}
			else {
				break;
			}
		}
	}
	totalCandidatePatterns += candiColocCounter;
}

void colocationFinder::candidateColocationGeneral(Integer degree) {
	Integer lastIndex = degree - 2;
	bool flag;
	candiColocCounter = 0;
	h_candiColocations.clear();
	h_candiColocations.shrink_to_fit();
	for (size_t i = 0; i < h_prevelentColocationCount; i++) {
		flag = true;
		for (size_t j = i + 1; j < h_prevelentColocationCount; j++) {
			for (size_t k = 0; k < lastIndex; k++) {
				if (h_prevalantColocations[i*(degree - 1) + k] != h_prevalantColocations[j*(degree - 1) + k]) {
					flag = false;
					break;
				}
			}
			if (flag && h_prevalantColocations[i*(degree - 1) + lastIndex]< h_prevalantColocations[j*(degree - 1) + lastIndex]) {
				vector<Integer> inter;
				for (Integer t = 0; t < degree - 1; t++) {
					inter.push_back(h_prevalantColocations[i*(degree - 1) + t]);
				}
				inter.push_back(h_prevalantColocations[j*(degree - 1) + lastIndex]);
				//module to generate k-1 degree subset
				flag = generateandCheckSubsets(inter, degree);
				clearSubsetVectors();
				if (!flag) {
					break;
				}
				else {
					for (size_t l = 0; l < degree; l++) {
						h_candiColocations.push_back(inter[l]);
					}
					candiColocCounter++;
				}
			}
			else {
				break;
			}
		}
	}
	totalCandidatePatterns += candiColocCounter;
}
void colocationFinder::kplus2ColocationGeneral(Integer degree) {
	Integer lastIndex = degree - 2;
	bool flag;
	candiColocCounter = 0;
	kplus2CandiColocation.clear();
	kplus2CandiColocation.shrink_to_fit();
	for (size_t i = 0; i < h_prevelentColocationCount2; i++) {
		flag = true;
		for (size_t j = i + 1; j < h_prevelentColocationCount2; j++) {
			for (size_t k = 0; k < lastIndex; k++) {
				if (h_prevalantColocations2[i*(degree - 1) + k] != h_prevalantColocations2[j*(degree - 1) + k]) {
					flag = false;
					break;
				}
			}
			if (flag && h_prevalantColocations2[i*(degree - 1) + lastIndex]< h_prevalantColocations2[j*(degree - 1) + lastIndex]) {
				vector<Integer> inter;
				for (Integer t = 0; t < degree - 1; t++) {
					inter.push_back(h_prevalantColocations2[i*(degree - 1) + t]);
				}
				inter.push_back(h_prevalantColocations2[j*(degree - 1) + lastIndex]);
				//module to generate k-1 degree subset
				flag = generateandCheckSubsetskplus2(inter, degree);
				clearSubsetVectors();
				if (!flag) {
					break;
				}
				else {
					for (size_t l = 0; l < degree; l++) {
						kplus2CandiColocation.push_back(inter[l]);
					}
				}
			}
			else {
				break;
			}
		}
	}
}

void colocationFinder::subsetGen(vector<Integer> &inter, Integer k, Integer n, Integer idx) {
	if (idx == n)
		return;

	if (k == 1) {
		for (size_t i = idx; i<n; i++)
		{
			subset.push_back(inter[i]);
			subsetList.push_back(subset);
			subset.pop_back();
		}
	}

	for (size_t j = idx; j<n; j++) {
		subset.push_back(inter[j]);
		subsetGen(inter, k - 1, n, j + 1);
		subset.pop_back();
	}
}

bool colocationFinder::checkSubset(vector<Integer> subsetElem) {
	Integer degree = subsetElem.size();
	Integer flag = true;
	for (size_t i = 0; i < h_prevelentColocationCount; i++) {
		for (size_t j = 0; j < degree; j++) {
			if (h_prevalantColocations[i*degree + j] != subsetElem[j]) {
				flag = false;
				break;
			}
			flag = true;
		}
		if (flag) {
			return flag;
		}
	}
	return flag;
}
bool colocationFinder::checkSubset2(vector<Integer> subsetElem) {
	Integer degree = subsetElem.size();
	Integer flag = true;
	for (size_t i = 0; i < candiColocCounter; i++) {
		for (size_t j = 0; j < degree; j++) {
			if (h_candiColocations[i*degree + j] != subsetElem[j]) {
				flag = false;
				break;
			}
			flag = true;
		}
		if (flag) {
			return flag;
		}
	}
	return flag;
}
bool colocationFinder::checkSubsetkplus2(vector<Integer> subsetElem) {
	Integer degree = subsetElem.size();
	Integer flag = true;
	for (size_t i = 0; i < h_prevelentColocationCount2; i++) {
		for (size_t j = 0; j < degree; j++) {
			if (h_prevalantColocations2[i*degree + j] != subsetElem[j]) {
				flag = false;
				break;
			}
			flag = true;
		}
		if (flag) {
			return flag;
		}
	}
	return flag;
}
bool colocationFinder::generateandCheckSubsets(vector<Integer> &inter, Integer degree) {
	subsetGen(inter, degree - 1, degree, 0);
	bool flag = true;
	for (size_t i = 2; i < degree; i++) {
		flag = checkSubset(subsetList[i]);
		if (!flag) {
			return flag;
		}
	}
	return flag;
}
bool colocationFinder::generateandCheckSubsets2(vector<Integer> &inter, Integer degree) {
	subsetGen(inter, degree - 1, degree, 0);
	bool flag = true;
	for (size_t i = 2; i < degree; i++) {
		flag = checkSubset2(subsetList[i]);
		if (!flag) {
			return flag;
		}
	}
	return flag;
}
bool colocationFinder::generateandCheckSubsetskplus2(vector<Integer> &inter, Integer degree) {
	subsetGen(inter, degree - 1, degree, 0);
	bool flag = true;
	for (size_t i = 2; i < degree; i++) {
		flag = checkSubsetkplus2(subsetList[i]);
		if (!flag) {
			return flag;
		}
	}
	return flag;
}

Integer colocationFinder::getIndex(vector<Integer> inner, Integer degree) {
	bool flag = false;
	for (size_t i = 0; i < h_prevelentColocationCount; i++) {
		for (size_t j = 0; j < degree; j++) {
			if (inner[j] != h_prevalantColocations[i*degree + j]) {
				flag = false;
				break;
			}
			flag = true;
		}
		if (flag) {
			return i;
		}
	}
	return -1;
}
Integer colocationFinder::getIndex2(vector<Integer> inner, Integer degree) {
	bool flag = false;
	for (size_t i = 0; i < h_prevelentColocationCount2; i++) {
		for (size_t j = 0; j < degree; j++) {
			if (inner[j] != h_prevalantColocations2[i*degree + j]) {
				flag = false;
				break;
			}
			flag = true;
		}
		if (flag) {
			return i;
		}
	}
	return -1;
}

void colocationFinder::generatePrevalentPatternsGeneral(Integer degree) {
	memoryTracker = 0.0;
	Integer lastIndex = degree - 2;
	h_prevelentColocationCount2 = 0;
	auto start = get_time::now();
	Integer filterpruneCounter = 0;
	if (estimatePrevalent < candiColocCounter) {
		free(h_prevalentSlotCounts);
		h_prevalentSlotCounts = (Integer**)malloc((candiColocCounter)*sizeof(Integer*));
		estimatePrevalent = candiColocCounter;
		for (int i = 0; i < estimatePrevalent; i++) {
			h_prevalentSlotCounts[i] = (Integer*)malloc(slotLimit*sizeof(Integer));
		}
	}
	for (size_t i = 0; i < candiColocCounter; i++) {
		for (Integer t = 0; t < degree; t++) {
			h_coloc[t] = h_candiColocations[i*degree + t];
		}
		//Filtering Mechanism
		if (parameter.FilterON_OFF == 1) {
			auto fstart = get_time::now();
			Real upperbound = filterCPP(h_coloc, degree);
			if (upperbound < parameter.PIthreshold) {
				filterpruneCounter++;
				//std::cout << "Degree " << degree << ": Candidate pattern " << i + 1 << "/" << candiColocCounter << " end (filter based pruned)" << std::endl;
				auto fend = get_time::now();
				auto fdiff = fend - fstart;
				Real msecs = chrono::duration_cast<chrono::milliseconds>(fdiff).count();
				totalFilterTime += (msecs / 1000.0);
				continue;
			}
			auto fend = get_time::now();
			auto fdiff = fend - fstart;
			Real msecs = chrono::duration_cast<chrono::milliseconds>(fdiff).count();
			totalFilterTime += (msecs / 1000.0);
		}
		//Filtering Mechanism Ends
		auto rStart = get_time::now();
		Integer table1Index;
		Integer table2Index;
		vector<Integer> combiColoc1;    //k combining pattern to make a k+1 candidate pattern 
		vector<Integer> combiColoc2;
		Integer instaceCountTable1;
		for (size_t t = 0; t < degree; t++) {
			if (t < degree - 2) {
				combiColoc1.push_back(h_candiColocations[i*degree + t]);
				combiColoc2.push_back(h_candiColocations[i*degree + t]);
			}
		}
		combiColoc1.push_back(h_candiColocations[i*degree + (degree - 2)]);
		combiColoc2.push_back(h_candiColocations[i*degree + (degree - 1)]);

		table1Index = getIndex(combiColoc1, degree - 1);
		instaceCountTable1 = h_prevalentInstanceTableSize[table1Index];///(degree-1);
		//std::cout << "Table1 size =" << instaceCountTable1 << "h_slot+indexes :" << 2 * instaceCountTable1 << std::endl;
		if (instaceCountTable1 > slotLimit) {
			free(h_prevalentSlotCounts[h_prevelentColocationCount2]);
			h_prevalentSlotCounts[h_prevelentColocationCount2] = (Integer*)malloc(instaceCountTable1*sizeof(Integer));
			slotLimit = instaceCountTable1;
		}
		memset(h_prevalentSlotCounts[h_prevelentColocationCount2], 0, sizeof(Integer)*instaceCountTable1);
		Integer bitmapSize = 0;
		Integer fType;
		for (size_t Idx = 0; Idx < degree; Idx++) {
			fType = h_coloc[Idx];
			bitmapSize += featureInstanceCount[fType];
		}
		vector<Integer> bitmap(bitmapSize, 0);
		Integer lastFeatureID = h_coloc[degree - 1];
		generalInstanceTableSlotCounter_Kernel(instaceCountTable1, h_prevalentSlotCounts[h_prevelentColocationCount2], degree, h_prevalentInstanceTable[table1Index], bitmap, h_coloc, lastFeatureID);

		totalKernelLaunches++;
			Real PI = getPI(bitmap, degree, h_coloc);
			if (PI < parameter.PIthreshold) {
			auto rend = get_time::now();
			auto rdiff = rend - rStart;
			Real mrsecs = chrono::duration_cast<chrono::milliseconds>(rdiff).count();
			totalRefineTime += (mrsecs / 1000.0);
			continue;
		}
		for (Integer Idx = 0; Idx < degree; Idx++) {
			h_prevalantColocations2.push_back(h_coloc[Idx]);
		}
		h_prevelentColocationCount2++;
		auto rend = get_time::now();
		auto rdiff = rend - rStart;
		Real mrsecs = chrono::duration_cast<chrono::milliseconds>(rdiff).count();
		totalRefineTime += (mrsecs / 1000.0);
	}
	totalFilteredPatterns += filterpruneCounter;
}

void colocationFinder::generateInstanceTableGeneral(Integer degree) {
	//std::cout << "Degree = " << degree<<std::endl;
	Integer lastIndex = degree - 2;
	h_prevalentInstanceTable2 = (Integer**)malloc((h_prevelentColocationCount2)*sizeof(Integer*));
	StatusType status;
	auto rStart = get_time::now();
	for (size_t i = 0; i < h_prevelentColocationCount2; i++) {
		if (hasInstanceTable[i]==1)
		{
			Integer instaceCountTable1;
			Integer table1Index;
			Integer table2Index;
			vector<Integer> combiColoc1;    //k combining pattern to make a k+1 candidate pattern 
			vector<Integer> combiColoc2;

			for (size_t t = 0; t < degree; t++) {
				if (t < degree - 2) {
					combiColoc1.push_back(h_prevalantColocations2[i*degree + t]);
					combiColoc2.push_back(h_prevalantColocations2[i*degree + t]);
				}
			}
			combiColoc1.push_back(h_prevalantColocations2[i*degree + (degree - 2)]);
			combiColoc2.push_back(h_prevalantColocations2[i*degree + (degree - 1)]);

			table1Index = getIndex(combiColoc1, degree - 1);
			instaceCountTable1 = h_prevalentInstanceTableSize[table1Index];
			//std::cout << "instaceCountTable1 = " << instaceCountTable1 << std::endl;
			size_t totalInstances = 0;
			if (instaceCountTable1 > instanceMax) {
				free(h_Indexes);
				h_Indexes = (Index64*)malloc(instanceMax*sizeof(Index64));
				}
			exclusive_scan(i, instaceCountTable1);
			totalInstances = getTotalInstances(h_prevalentSlotCounts[i], instaceCountTable1);
			size_t extraMemoryRequired = totalInstances * degree * 4;
			Real mem = (Real)(extraMemoryRequired) / (1024 * 1024 * 1024);
			memoryTracker += mem;
			//std::cout << "degree =" << degree << " mem = " << mem << std::endl;
			#pragma omp parallel for
			for (int TID = 0; TID < instaceCountTable1; TID++) {
				size_t temp = h_Indexes[TID];
				size_t value = temp * degree;
				h_Indexes[TID] = value;
			}
			Integer lastFeatureID = h_prevalantColocations2[i*degree + (degree - 1)];
				h_prevalentInstanceTableSize2.push_back(totalInstances);
				size_t mSize = totalInstances*degree;
				h_prevalentInstanceTable2[i] = (Integer*)malloc(mSize*sizeof(Integer));
				generateInstanceTableGeneral_Kernel(instaceCountTable1, h_prevalentSlotCounts[i], h_Indexes, h_prevalentInstanceTable2[i], degree, h_prevalentInstanceTable[table1Index], 0, lastFeatureID);
				totalKernelLaunches++;
		}
		else {
			h_prevalentInstanceTableSize2.push_back(0);
		}
	}
	auto rend = get_time::now();
	auto rdiff = rend - rStart;
	Real mrsecs = chrono::duration_cast<chrono::milliseconds>(rdiff).count();
	totalRefineTime += (mrsecs / 1000.0);
}

void colocationFinder::copyPrevalentColocations() {
	h_prevalantColocations = h_prevalantColocations2;
	h_prevalantColocations2.clear();
	h_prevalantColocations2.shrink_to_fit();
	h_prevalentInstanceTable = h_prevalentInstanceTable2;
	h_prevelentColocationCount = h_prevelentColocationCount2;
	h_prevelentColocationCount2 = 0;
	h_prevalentInstanceTableSize = h_prevalentInstanceTableSize2;
	h_prevalentInstanceTableSize2.clear();
	h_prevalentInstanceTableSize2.shrink_to_fit();
}
//Cleaning Functions start
void colocationFinder::resetPrevalentData(Integer degree) {

}

void colocationFinder::clearMemory() {
	free(h_prevalentInstanceTable);

	h_prevalantColocations.clear();
	h_prevalantColocations.shrink_to_fit();

	h_prevalantColocations2.clear();
	h_prevalantColocations2.shrink_to_fit();

	h_prevalentInstanceTableSize.clear();
	h_prevalentInstanceTableSize.shrink_to_fit();

	h_prevalentInstanceTableSize2.clear();
	h_prevalentInstanceTableSize2.shrink_to_fit();
	delete[] h_coloc;
}

void colocationFinder::clearSubsetVectors() {
	subset.clear();
	subset.shrink_to_fit();
	subsetList.clear();
	subsetList.shrink_to_fit();
}

//Cleaning Functions ends

//Log Function starts
void colocationFinder::savetoFile2(Integer degree) {
	std::string degStr = std::to_string(degree);
	std::string FileName = location + "ColocationDegree_" + degStr + ".txt";
	std::ofstream fout(FileName.c_str());
	for (size_t i = 0; i < h_prevelentColocationCount; i++) {
		fout << "( ";
		for (size_t j = 0; j < degree; j++) {
			size_t index = h_prevalantColocations[i*degree + j];
			fout << featureTypes[index] << " ";
			if (j != degree - 1) {
				fout << "| ";
			}
			else {
				fout << ")" << "\n";
			}
		}
	}
	fout.close();
}
void colocationFinder::savetoFileGen(Integer degree) {
	std::string degStr = std::to_string(degree);
	std::string FileName = location + "ColocationDegree_" + degStr + ".txt";
	std::ofstream fout(FileName.c_str());
	for (size_t i = 0; i < h_prevelentColocationCount2; i++) {
		fout << "( ";
		for (size_t j = 0; j < degree; j++) {
			size_t index = h_prevalantColocations2[i*degree + j];
			fout << featureTypes[index] << " ";
			if (j != degree - 1) {
				fout << "| ";
			}
			else {
				fout << ")" << "\n";
			}
		}
	}
	fout.close();
}

void colocationFinder::timeLog(Integer degree, std::string function, std::string eventType, Integer duration) {
	time_t currentTime;
	struct tm *localTime;
	time(&currentTime);                   // Get the current time
	localTime = localtime(&currentTime);  // Convert the current time to the local time
	timeLogStream.open(logFileName.c_str(), ofstream::app);
	timeLogStream << degree << " " << function << " " << eventType << " " << localTime->tm_year + 1900 << "\\" << localTime->tm_mon + 1 << "\\" << localTime->tm_mday << " " << localTime->tm_hour << ":" << localTime->tm_min << ":" << localTime->tm_sec << "\n";
	if (eventType == "end") {
		timeLogStream << degree << " " << function << " " << "Duration: " << " " << duration << "seconds \n";
	}
	timeLogStream.close();
}
void colocationFinder::compLog(Integer degree, std::string function) {
	compLogStream.open(compLogFileName.c_str(), ofstream::app);
	compLogStream << "Degree: " << degree << ", " << function << "\n";
	compLogStream.close();
}
void colocationFinder::compLog(Integer i, Integer degree, Integer candicolocNumber, std::string Remark, std::string PIType, Real PI) {
	compLogStream.open(compLogFileName.c_str(), ofstream::app);
	if (Remark != "") {
		compLogStream << "Degree: " << degree << ", Candidate " << i << "/" << candicolocNumber << PIType << PI << " (" << Remark << ") \n";
	}
	else {
		compLogStream << "Degree: " << degree << ", Candidate " << i << "/" << candicolocNumber << PIType << PI << "\n";
	}
	compLogStream.close();
}
void colocationFinder::compLog(Integer degree, Integer i, std::string totalmem, Integer candiColocNum) {
	compLogStream.open(compLogFileName.c_str(), ofstream::app);
	compLogStream << "Degree: " << degree << " pattern No.: " << i << "/" << candiColocNum << " Total memRequired =" << totalmem << "GB\n";
	compLogStream.close();
}
void colocationFinder::loggingFileOpenFunction() {
	logFileName = location + "GPU_TimeLog.txt";
	compLogFileName = location + "GPU_ComputationLog.txt";
}
void colocationFinder::loggingFileCloseFunction() {
	timeLogStream.open(logFileName.c_str(), ofstream::app);
	timeLogStream << "distance Threshold: " << parameter.thresholdDistance << "\n";
	timeLogStream << "prevalance Threshold: " << parameter.PIthreshold << "\n";
	timeLogStream.close();
}
//Log Functions ends
void colocationFinder::calcInstanceTablePatternIndexes(Integer degree) {
	hasInstanceTable.clear();
	hasInstanceTable.shrink_to_fit();
	for (Integer i = 0; i < h_prevelentColocationCount2; i++) {
		hasInstanceTable.push_back(0);
	}
	Integer degreePlus1 = degree + 1;
	kplus2ColocationGeneral(degree+1);
	Integer kplu2PatternsCount = kplus2CandiColocation.size()/ degreePlus1;
	if (kplu2PatternsCount == 0) {
		needInstanceTable = false;
		return;
	}
	for (Integer i = 0; i < kplu2PatternsCount; i++) {
		vector<Integer> combiColoc1;
		for (size_t t = 0; t < degreePlus1; t++) {
			if (t < degreePlus1 - 2) {
				combiColoc1.push_back(kplus2CandiColocation[i*degreePlus1 + t]);
				}
		}
		combiColoc1.push_back(kplus2CandiColocation[i*degreePlus1 + (degreePlus1 - 2)]);

		Integer table1Index = getIndex2(combiColoc1, degreePlus1 - 1);
		hasInstanceTable[table1Index] = 1;
	}
}
void colocationFinder::Begin(Integer argc, char**argv) {
	auto begin = get_time::now();
	std::string datasetFilename;
	if (argc > 1) {
		std::string configurationFile = argv[1];
		ifstream in(configurationFile.c_str());
		std::string line;
		std::string parameterFilename;
		getline(in, line);
		location = line;
		getline(in, line);
		parameterFilename = line;
		getline(in, line);
		datasetFilename = line;
		getline(in, line);
		outputLocation = line;
		getline(in, line);
		outputFile = line;
		in.close();
		LoadParameter(parameterFilename);
		auto start = get_time::now();
		//timeLog(1, "populate Data", "start", 0);
		populateData(datasetFilename);
		auto end = get_time::now();
		auto diff = end - start;
		//timeLog(1, "populate Data", "end", chrono::duration_cast<sec>(diff).count());
		//std::cout << "Data Load time (PopulateData function) :  " << chrono::duration_cast<sec>(diff).count() << " seconds " << endl;

	}
	else {
		//std::cout << "Configuration file missing...";
		exit(0);
	}
	auto processingStart = get_time::now();
	degree2Processing();
	std::cout << "Saving degree 2 colcation data..." << std::endl;
	savetoFile2(2);
	std::cout << "\n\nDegree 2 processing ends...\n" << std::endl;
	Integer degree = 3;
	Integer featureCount = featureTypes.size();
	tableGenRequired(degree+1);
	while (candiColocCounter > 0 && degree <= featureCount) {
		std::cout << "\n\nDegree " << degree << " processing starts..." << std::endl;
		//general Table Instance Generation
		generatePrevalentPatternsGeneral(degree);
		if (needInstanceTable) {
			calcInstanceTablePatternIndexes(degree);
			generateInstanceTableGeneral(degree);
		}
		//std::cout << "Saving degree " << degree << " colcation data..." << std::endl;
		savetoFileGen(degree);
		//std::cout << "\n\nDegree " << degree << " processing ends...\n" << std::endl;
		copyPrevalentColocations();

		//candidate colocation pattern generator 
		//std::cout << "\n\nDegree " << degree + 1 << " processing starts..." << std::endl;
		auto start = get_time::now();
		//timeLog(degree + 1, "candidateColocationGeneral ", "start", 0);
		candidateColocationGeneral(degree + 1);
		tableGenRequired(degree+2);
		auto end = get_time::now();
		auto diff = end - start;
		//timeLog(degree + 1, "candidateColocationGeneral " + std::to_string(candiColocCounter), "end", chrono::duration_cast<sec>(diff).count());
		//std::cout << "candidateColocationGeneral() took  " << chrono::duration_cast<sec>(diff).count() << " seconds " << endl;

		std::cout << "\n\nDegree " << degree << " processing ends..." << std::endl;
		if (candiColocCounter == 0) {
			clearMemory();//clean memory
			break;
		}
	
		degree++;
	}
	auto last = get_time::now();
	auto totalTimeDifference = last - begin;  //ioTime + preprocessing time + processing time
	auto pt = last - processingStart;
	Real processingTime = chrono::duration_cast<chrono::milliseconds>(pt).count();
	std::cout << std::endl << "Processing Time : " << processingTime << " milliseconds " << endl;
	Real totalTime = chrono::duration_cast<chrono::milliseconds>(totalTimeDifference).count();
	std::ofstream expResult;
	std::string expFile = outputLocation + outputFile;
	expResult.open(expFile.c_str(), ofstream::app);
	std::string exe = argv[0];
	int pos = 0;
	for (int i = 0; i < exe.length(); i++) {
		if (exe[i] == '\\') {
			pos = i;
		}
	}
	exe.erase(0, pos + 1);
	int pos2 = 0;
	for (int i = 0; i < datasetFilename.length(); i++) {
		if (datasetFilename[i] == '\\') {
			pos2 = i;
		}
	}
	datasetFilename.erase(0, pos2 + 1);
	expResult << datasetFilename << "," << exe << "," << parameter.PIthreshold << "," << "1" << "," << ioTime << "," << preprocessingTime << "," << totalFilterTime << "," << preprocessingTime- totalFilterTime << "," << 0 << "," << totalTime  << "," << processingTime+preprocessingTime << "\n";

	expResult.close();

}





Real colocationFinder::filterCPP(Integer * h_coloc, Integer degree) {
	//Integer *countMap;
	vector<Integer> countMap(degree* gridStructure.totalCells, 0);
	Integer tempsum = 0;
	#pragma omp parallel for
	for (int TID = 0; TID < gridStructure.totalCells; TID++) {
		Integer neighborNumber;
		Integer neighborcells[4];
		neighborNumber = 0;
		Integer Range = 1;
		if ((TID + gridStructure.gridSize.x) < gridStructure.totalCells) {
			Integer CIDY = TID / gridStructure.gridSize.x;
			Integer CIDX = TID - (gridStructure.gridSize.x*CIDY);
			for (Integer j = 0; j <= Range; ++j)
			{
				for (Integer i = 0; i <= Range; ++i)
				{
					SInteger NCID = getNeighborCellID(CIDX + i, CIDY + j, TID);
					if (NCID < 0 || NCID >= gridStructure.totalCells)
					{
						continue;
					}
					neighborcells[neighborNumber] = NCID;
					neighborNumber++;
				}
			}
			if (neighborNumber < 4) {
				continue;
			}
			checkandUpdate(degree, h_coloc, neighborcells, countMap);
		}
	}

	Real PI = 1.0;
	Integer from;
	Integer to;
	for (Integer j = 0; j < degree; j++) {
		Integer index = h_coloc[j];
		Integer totalInstance = featureInstanceCount[index];
		from = j*gridStructure.totalCells;
		to = (j + 1)*gridStructure.totalCells - 1;
		size_t bitSum = std::accumulate(countMap.begin() + from, countMap.begin() + to, 0);
			Real pr = bitSum / (Real)totalInstance;
		if (pr < PI) {
			PI = pr;
		}
	}
	countMap.clear();
	countMap.shrink_to_fit();
	return PI;
}

void colocationFinder::degree2TableSlotCounter_Kernel(size_t maxCount, Integer* slots, size_t start, size_t end, size_t secondstart, vector<Integer>& d_bitmap, Integer combiningFeatureID)
{
#pragma omp parallel for
	for (int TID = 0; TID < maxCount; TID++) {
		size_t idx1 = TID + start;
		Integer currentCell = instanceCellIDs[idx1];
		Integer CIDY = currentCell / gridStructure.gridSize.x;
		Integer CIDX = currentCell - (gridStructure.gridSize.x*CIDY);
		SInteger Range = 1;
		for (SInteger j = -1; j <= Range; ++j)
		{
			for (SInteger i = -1; i <= Range; ++i)
			{
				SInteger NCID = getNeighborCellID(CIDX + i, CIDY + j, currentCell);
				if (NCID < 0 || NCID >= gridStructure.totalCells)
				{
					continue;
				}
				size_t index = NCID*maxFeaturesNum + combiningFeatureID;
				size_t count = cellEventInstanceCount[index];
				if (count > 0) {
					size_t startIndex = cellEventInstanceStart[index];
					size_t endIndex = cellEventInstanceEnd[index];
					for (size_t pos = startIndex; pos <= endIndex; pos++) {
						Integer instanceID = cellBasedSortedIndex[pos];
						Real dist = DistanceInMeters(idx1, instanceID);
						if (dist <= parameter.squaredDistanceThreshold) {
							Integer slotValue = slots[TID];
							slotValue++;
							slots[TID] = slotValue;
							size_t firstcount = end - start;
							size_t bitIndex1 = TID;
							size_t bitIndex2 = instanceID - secondstart + firstcount;
							d_bitmap[bitIndex1] = 1;
							d_bitmap[bitIndex2] = 1;
						}
					}
				}
			}
		}
	}
}

void colocationFinder::scalebyConstant_Kernel(Integer degree, Index64* indexes, size_t instaceCountTable1) {
	for (size_t TID = 0; TID < instaceCountTable1; TID++) {
		size_t temp = indexes[TID];
		size_t value = temp * degree;
		indexes[TID] = value;
	}
}

void colocationFinder::degree2TableGenerator_Kernel(Integer table1InstanceCount, Integer*slots, Index64* indexes, Integer*& d_intermediate, size_t start, size_t end, Integer combiningFeatureID)
{
	//size_t idx1 = 0;
#pragma omp parallel for
	for (int TID = 0; TID  < table1InstanceCount; TID++) {
		size_t idx1 = TID + start;
		Integer slotsValue = slots[TID];
		if (slotsValue == 0) {
			continue;
		}
		size_t index = indexes[TID];
		Integer currentCell = instanceCellIDs[idx1];
		Integer CIDY = currentCell / gridStructure.gridSize.x;
		Integer CIDX = currentCell - (gridStructure.gridSize.x*CIDY);

		SInteger Range = 1;
		for (SInteger j = -1; j <= Range; ++j)
		{
			for (SInteger i = -1; i <= Range; ++i)
			{
				SInteger NCID = getNeighborCellID(CIDX + i, CIDY + j, currentCell);
				if (NCID < 0 || NCID >= gridStructure.totalCells)
				{
					continue;
				}
				size_t indx = NCID*maxFeaturesNum + combiningFeatureID;
				size_t count = cellEventInstanceCount[indx];
				if (count > 0) {
					size_t startIndex = cellEventInstanceStart[indx];
					size_t endIndex = cellEventInstanceEnd[indx];
					for (size_t pos = startIndex; pos <= endIndex; pos++) {
						Integer instanceID = cellBasedSortedIndex[pos];
						Real dist = DistanceInMeters(idx1, instanceID);
						if (dist <= parameter.squaredDistanceThreshold) {
							d_intermediate[index] = idx1;
							index++;
							d_intermediate[index] = instanceID;
							index++;
						}
					}
				}
			}
		}
	}
}


void colocationFinder::checkandUpdate(Integer degree, Integer* d_coloc, Integer* quadrant, vector<Integer>& countMap) {
	Integer check = 0;
	for (Integer i = 0; i < degree; i++) {
		Integer eventID = d_coloc[i];
		for (Integer j = 0; j < 4; j++) {
			Integer cellID = quadrant[j];
			Integer index = cellID*maxFeaturesNum + eventID;
			if (cellEventInstanceCount[index]>0) {
				check++;
				break;
			}
		}
		if (check <= i) {
			break;
		}
	}
	if (check != degree) {
		return;
	}
	for (Integer i = 0; i < degree; i++) {
		for (Integer j = 0; j < 4; j++) {
			Integer cellID = quadrant[j];
			Integer eventID = d_coloc[i];
			Integer index = cellID*maxFeaturesNum + eventID;
			if (cellEventInstanceCount[index]>0) {
				Integer value = cellEventInstanceCount[index];
				countMap[i*gridStructure.totalCells + cellID] = value;
			}
		}
	}
}

Real colocationFinder::getPI(vector<Integer> bitmap, Integer degree, Integer *coloc) {
	Real PI = 1.0;
	Integer from;
	Integer to;
	for (Integer j = 0; j < degree; j++) {
		Integer index = coloc[j];
		Integer totalInstance = featureInstanceCount[index];
		if (j == 0) {
			from = 0;
			to = totalInstance;
		}
		else {
			from = to - 1;
			to = from + totalInstance + 1;
		}
		size_t bitSum = std::accumulate(bitmap.begin() + from, bitmap.begin() + to - 1, 0);
		Real pr = bitSum / (Real)totalInstance;
		if (pr < PI) {
			PI = pr;
		}
	}
	return PI;

}


void colocationFinder::generalInstanceTableSlotCounter_Kernel(size_t instaceCountTable1, Integer* slots, Integer degree, Integer* table1, vector<Integer>& d_bitmap, Integer* coloc, Integer lastFeatureID)
{
	#pragma omp parallel for
	for (int TID = 0; TID < instaceCountTable1; TID++) {
		Integer degreek = degree - 1;
		size_t counter = 0;
		SInteger Range = 1;
		Integer currentInstanceID = table1[TID*degreek];
		Integer currentCell = instanceCellIDs[currentInstanceID];
		Integer CIDY = currentCell / gridStructure.gridSize.x;
		Integer CIDX = currentCell - (gridStructure.gridSize.x*CIDY);
		for (SInteger j = -1; j <= 1; ++j) {
			for (SInteger i = -1; i <= 1; ++i) {
				SInteger NCID = getNeighborCellID(CIDX + i, CIDY + j, currentCell);
				if (NCID < 0 || NCID >= gridStructure.totalCells) {
					continue;
				}
				size_t lastIndex = NCID*maxFeaturesNum + lastFeatureID;
				size_t count = cellEventInstanceCount[lastIndex];
				if (count > 0) {
					size_t startIndex = cellEventInstanceStart[lastIndex];
					size_t endIndex = cellEventInstanceEnd[lastIndex];
					for (size_t pos = startIndex; pos <= endIndex; pos++) {
						Integer lastInstanceID = cellBasedSortedIndex[pos];
						Integer flag2 = 1;
						for (int degIndex = 0; degIndex < degreek; degIndex++) {
							Integer cIndex = TID*degreek + degIndex;
							currentInstanceID = table1[cIndex];
							Real dist = DistanceInMeters(currentInstanceID, lastInstanceID);
							if (dist > parameter.squaredDistanceThreshold) {
								flag2 = 0;
								break;
							}
						}
						if (flag2 == 1) {
							counter++;
							slots[TID] = counter;
							size_t offset = 0;
							Integer bitIndex = 0;
							Integer fType;
							Integer value = 0;
							for (size_t idx4 = 0; idx4 < degreek; idx4++) {
								size_t idxCurrent = TID*degreek + idx4;
								value = table1[idxCurrent];
								fType = coloc[idx4];
								size_t fstart = featureInstanceStart[fType];
								size_t fcount = featureInstanceCount[fType];
								bitIndex = value - fstart + offset;
								offset += fcount;
								d_bitmap[bitIndex] = 1;
							}
							Integer lastStart = featureInstanceStart[lastFeatureID];
							bitIndex = lastInstanceID - lastStart + offset;
							d_bitmap[bitIndex] = 1;
						}
					}
				}
			}
		}
	}
}


void colocationFinder::generateInstanceTableGeneral_Kernel(size_t instaceCountTable1, Integer* slots, Index64* indexes, Integer* d_intermediate, Integer degree, Integer* table1, size_t startFrom, Integer lastFeatureID)
{
	#pragma omp parallel for
	for (int TID = 0; TID < instaceCountTable1; TID++) {
		Integer slotsValue = slots[TID];
		if (slotsValue == 0) {
			continue;
		}

		size_t index;
		if (startFrom > 0) {
			size_t indexPos = indexes[startFrom];
			index = indexes[TID] - indexPos;
		}
		else {
			index = indexes[TID];
		}
		Integer degreek = degree - 1;
		SInteger Range = 1;
		Integer currentInstanceID = table1[TID*degreek];
		Integer currentCell = instanceCellIDs[currentInstanceID];
		Integer CIDY = currentCell / gridStructure.gridSize.x;
		Integer CIDX = currentCell - (gridStructure.gridSize.x*CIDY);
		for (SInteger j = -1; j <= Range; ++j)
		{
			for (SInteger i = -1; i <= Range; ++i)
			{
				SInteger NCID = getNeighborCellID(CIDX + i, CIDY + j, currentCell);
				if (NCID < 0 || NCID >= gridStructure.totalCells)
				{
					continue;
				}
				size_t lastFeatureCellIndex = NCID*maxFeaturesNum + lastFeatureID;
				size_t count = cellEventInstanceCount[lastFeatureCellIndex];
				if (count > 0) {
					size_t startIndex = cellEventInstanceStart[lastFeatureCellIndex];
					size_t endIndex = cellEventInstanceEnd[lastFeatureCellIndex];
					for (size_t pos = startIndex; pos <= endIndex; pos++) {
						Integer lastInstanceID = cellBasedSortedIndex[pos];
						Integer flag2 = 1;
						for (int degIndex = 0; degIndex < degreek; degIndex++) {
							currentInstanceID = table1[TID*degreek + degIndex];
							Real dist = DistanceInMeters(currentInstanceID, lastInstanceID);
							if (dist > parameter.squaredDistanceThreshold) {
								flag2 = 0;
								break;
							}
						}
						if (flag2 == 1) {
							for (size_t idx4 = 0; idx4 < degreek; idx4++) {
								size_t idxCurrent = TID*degreek + idx4;
								Integer value = table1[idxCurrent];
								d_intermediate[index] = value;
								index++;
							}
							d_intermediate[index] = lastInstanceID;
							index++;
						}
					}
				}
			}
		}
	}
}

void colocationFinder::exclusive_scan(Integer index, Integer limit) {
	h_Indexes[0] = 0;
	for (int i = 0; i < limit - 1; i++) {
		h_Indexes[i + 1] = h_Indexes[i] + h_prevalentSlotCounts[index][i];
	}
}
Integer colocationFinder::getTotalInstances(Integer *slots, size_t table1InstanceCount) {
	Integer totalInstance = std::accumulate(slots, slots + table1InstanceCount, 0);
	return totalInstance;
}