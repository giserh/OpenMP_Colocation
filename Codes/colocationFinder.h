#pragma once
#include "DataTypes.h"
#include <limits>
#include <string>
#include <chrono>
#include <math.h>
#include <iomanip>
#include <ctime>
using  sec = chrono::seconds;
using get_time = chrono::steady_clock;
class colocationFinder {

public:
	colocationFinder(void);

public:
	~colocationFinder(void);

private:
	//Host Variables
	Real ioTime;
	Real preprocessingTime;
	size_t maxInstancesNum;
	size_t maxFeaturesNum;
	size_t maxCellNum;
	std::string location;
	vector <string> featureTypes;   //done
	Integer * instanceList;         //done
	Integer *featureInstanceStart;  //done
	Integer * featureInstanceEnd;	//done
	Integer * featureInstanceCount; //done
	Real * instanceLocationX;    //done
	Real * instanceLocationY;    //done
	Integer candiColocCounter;
	struct param parameter;
	//for subset
	vector<Integer> subset;
	vector<vector<Integer>> subsetList;
	vector<Integer> h_candiColocations;
	vector<Integer> kplus2CandiColocation;
	Integer *d_candiColocations;
	vector<Integer> tempColoc;

	Integer** h_prevalentInstanceTable;
	vector<Integer> h_prevalantColocations;
	Integer h_prevelentColocationCount;
	vector<Integer> h_prevalentInstanceTableSize;

	Integer** h_prevalentInstanceTable2;
	vector<Integer> h_prevalantColocations2;
	Integer h_prevelentColocationCount2;
	vector<Integer> h_prevalentInstanceTableSize2;
	vector<Integer> prevalanceCheck;
	Integer tempOldPrevalentColocCount;

	//added in version 2
	struct gridBox gridStructure;
	SInteger* cellEventInstanceCount;
	Integer * instanceCellIDs;

	//added in version 3
	SInteger* cellEventInstanceStart;
	SInteger* cellEventInstanceEnd;
	Integer* cellBasedSortedIndex;

	//for MRF
	vector<Integer> coarseInstanceStart;  //done
	vector<Integer> coarseInstanceEnd;	//done
	vector<Integer> coarseInstanceList;
	vector<vector<Integer>> coarseCandidateInstanceTables;
	vector<Integer> coarsePrevalance;
	Integer ** coarsePrevalentInstanceTables;
	Integer *coarsePrevalentInstanceTableSize;
	//Integer** cellEventInstanceCount;
	vector<Integer> coarseInstanceCount;

	//CUDA Device variables
	Integer *d_instanceList;
	Integer *d_featureInstanceStart;
	Integer *d_featureInstanceEnd;
	Integer *d_featureInstanceCount;
	Real	*d_instanceLocationX;
	Real	*d_instanceLocationY;
	SInteger* d_cellEventInstanceCount;

	SInteger* d_cellEventInstanceStart;
	SInteger* d_cellEventInstanceEnd;
	Integer* d_instanceCellIDs;
	Integer* d_cellBasedSortedIndex;

	//LOG VARIABLES
	std::string logFileName;
	std::ofstream timeLogStream;
	std::string compLogFileName;
	std::ofstream compLogStream;
	//std::string outputLocation;
	Real totalMemUsed;
	Real memoryTracker;
	std::string outputLocation;
	std::string outputFile;
	Integer totalCandidatePatterns;
	Integer totalFilteredPatterns;
	Real totalFilterTime;
	Real totalRefineTime;
	Real degree2FilterTime;

	bool needInstanceTable;
	Integer** h_prevalentSlotCounts;
	vector<Integer> hasInstanceTable;


	Integer *d_coloc;
	Integer *h_coloc;
	Integer *d_slotCounts;
	Index64 *h_Indexes;
	Index64 *d_Indexes;
	Index64 instanceMax;
	Integer *d_intermediate;
	Integer estimatePrevalent;
	Integer *d_table1;
	Integer slotLimit;
	Integer intermediateMax;
	Integer table1Max;
	Integer mapMax;
	Integer *d_countMap;
	//Integer *h_bitmap;
	Integer *d_bitmap;

	Real totalKernelTime;
	Integer totalKernelLaunches;
public:
	void Begin(Integer argc, char**argv);
	void candidateColocationGeneral(Integer degree);
	void populateData(std::string datasetFilename);
	void degree2CandidateColocationGenerator();
	Real DistanceInMeters(size_t idx1, size_t idx2);
	bool generateandCheckSubsets(vector<Integer> &inter, Integer degree);
	void subsetGen(vector<Integer> &inter, Integer k, Integer n, Integer idx);
	bool checkSubset(vector<Integer> subsetElem);
	void generateInstanceTableGeneral(Integer degree);
	void degree2Processing();
	void clearSubsetVectors();
	Integer getIndex(vector <Integer> inner, Integer degree);
	Integer getIndex2(vector <Integer> inner, Integer degree);
	void resetPrevalentData(Integer degree);
	void clearMemory();
	void copyPrevalentColocations();
	void LoadParameter(std::string parameterFileName);
	//void log(Integer degree, std::string function, std::string eventType);

	//new added in version 2
	void setGrid(Real maxX, Real maxY, Real minX, Real minY);
	Integer getCellID(Real lon, Real lat);
	//Integer getTypeID(size_t i);
	//bool checkifNeighbors(Integer cell1, Integer cell2);
	//void extractCoarsePrevalentColocations(Integer degree);
	
	
	//log related functions
	void timeLog(Integer degree, std::string function, std::string eventType, Integer duration);
	void loggingFileOpenFunction();
	void loggingFileCloseFunction();
	void compLog(Integer i, Integer degree, Integer candicolocNumber,std::string Remark, std::string PIType, Real PI);
	void compLog(Integer degree, std::string function);
	void compLog(Integer degree, Integer i, std::string totalmem, Integer candicolocNumber);
	void savetoFile2(Integer degree);
	void savetoFileGen(Integer degree);

	void generatePrevalentPatternsGeneral(Integer degree);
	void calcInstanceTablePatternIndexes(Integer degree);

	//test
	SInteger getNeighborCellID(Integer x, Integer y, Integer TID);
	void tableGenRequired(Integer degree);
	bool generateandCheckSubsets2(vector<Integer> &inter, Integer degree);
	bool checkSubset2(vector<Integer> subsetElem);
	void kplus2ColocationGeneral(Integer degree);
	bool checkSubsetkplus2(vector<Integer> subsetElem);
	bool generateandCheckSubsetskplus2(vector<Integer> &inter, Integer degree);

	Real filterCPP(Integer * h_coloc, Integer degree);
	Real getPI(vector<Integer> bitmap, Integer degree, Integer* coloc);
	void degree2TableSlotCounter_Kernel(size_t maxCount, Integer* slots, size_t start, size_t end, size_t secondstart, vector<Integer>& d_bitmap, Integer combiningFeatureID);
	void scalebyConstant_Kernel(Integer degree, Index64* indexes, size_t instaceCountTable1);
	void degree2TableGenerator_Kernel(Integer table1InstanceCount,Integer*slots, Index64* indexes, Integer*& d_intermediate, size_t start, size_t end, Integer combiningFeatureID);
	Real distanceinMeters(Integer idx1, Integer idx2);
	void checkandUpdate(Integer degree, Integer* d_coloc, Integer* quadrant, vector<Integer> & countMap);
	void generalInstanceTableSlotCounter_Kernel(size_t instaceCountTable1, Integer* slots, Integer degree, Integer* table1, vector<Integer>& d_bitmap, Integer* coloc, Integer lastFeatureID);
	void generateInstanceTableGeneral_Kernel(size_t instaceCountTable1, Integer* slots, Index64* indexes, Integer* d_intermediate, Integer degree, Integer* table1, size_t startFrom, Integer lastFeatureID);
	//void generalInstanceTableSlotCounter_Kernel(size_t instaceCountTable1, Integer* slots, Integer degree, Integer* table1, Integer* d_bitmap, Integer* coloc, Integer lastFeatureID)
	void exclusive_scan(Integer i, Integer limit);
	Integer getTotalInstances(Integer *slots, size_t table1InstanceCount);
}; 
