/*
 * megaMolSnapshotCommand.h
 *
 *  Created on: Mar 19, 2009
 *      Author: hpcdjenz
 */

#ifdef STEEREO
#ifndef MEGAMOLSNAPSHOTCOMMAND_H_
#define MEGAMOLSNAPSHOTCOMMAND_H_

#include <steereoCommand.h>
//#include <Simulation.h>

typedef struct {
	double boxLength;
	int numberOfComponents;
	int* numberOfParticles;
	//int* stepNumber;
} MS2DataType;

typedef struct {
	int npart;
	int npartMax;
} ComponentData;

class MegaMolSnapshotCommand : public SteereoCommand
{
	public:
  MegaMolSnapshotCommand ();
  ~MegaMolSnapshotCommand ();
  virtual ReturnType execute ();
  void setParameters (std::list<std::string> params);
  //static void setSimData (Simulation* simu) {sim = simu;};
  static SteereoCommand* generateNewInstance ();

  bool condition ();
  void setStepInterval (int interval) {stepInterval = interval;};


 private:
  // parameters needed for execution
  MS2DataType* lastInfo;
  static int startStep;
  int stepInterval;
  int colouringVal;

};

#endif /* MEGAMOLSNAPSHOTCOMMAND_H_ */
#endif /* STEEREO */
