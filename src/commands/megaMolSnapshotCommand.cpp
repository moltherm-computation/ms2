/*
 * megaMolSnapshotCommand.cpp
 *
 *  Created on: Mar 19, 2009
 *      Author: hpcdjenz
 */
#define STEEREO
#ifdef STEEREO

#include "megaMolSnapshotCommand.h"

#include <cmath>
#include <iostream>
#include <baseSimSteering.h>
#include <steereoStream.h>
#include <steereoDefinitions.h>
#include <steereoIntraCommunicator.h>


int MegaMolSnapshotCommand::startStep = -1;

MegaMolSnapshotCommand::MegaMolSnapshotCommand () : SteereoCommand (false, "getMegaMolSnapshot", 0)
{
  this->setCommandName ("getMegaMolSnapshot");
  stepInterval = 0;
  startStep = 0;
  lastInfo = NULL;
  SteereoLogger::setOutputLevel (2);
}

MegaMolSnapshotCommand::~MegaMolSnapshotCommand ()
{

}

ReturnType MegaMolSnapshotCommand::execute ()
{
  MS2DataType* ms2Data = (MS2DataType*) this->getData(0);
  logger->debug() << "The ms2Data pointer is " << ms2Data << std::endl;
  logger->debug() << "ms2Data->boxLength: " << ms2Data->boxLength << std::endl;
  logger->debug() << "ms2Data->numberOfComponents: " << ms2Data->numberOfComponents << std::endl;
  int numberOfComponents = (int) (ms2Data->numberOfComponents);
  int fieldSize, partNumber;

  SteereoStream processStream;
  SteereoStream* wholeStream;
  SteereoStream partStream[numberOfComponents];
  //SteereoStream outStream;
  //SteereoStream compStream[numberOfComponents];
  float minVal = 1e100, maxVal = 0;
  int counter = 0;

  if (this->getIntraCommunicator()->amIRoot())
  {
  	processStream.allocateMemory (8 * sizeof(float) + 3 * sizeof(int));
  	processStream << (float) -0.5 << (float) -0.5 << (float) -0.5;
  	processStream << (float) 0.5 << (float) 0.5 << (float) 0.5;
  	processStream << (int) this->getIntraCommunicator()->getSize();
  	processStream << (int) numberOfComponents;
  	processStream << (int) colouringVal;
    logger->debug() << "I have " << numberOfComponents << " components " << std::endl;
  }
  else
  {
  	processStream.allocateMemory (2 * sizeof(float));
  }

  for (int i = 0; i < numberOfComponents; i++)
  {
  	logger->debug() << "For component " << i << " I actually have " << ms2Data->numberOfParticles[i] << " particles" << std::endl;
  }


  for (int i = 0; i < numberOfComponents; i++)
  {
  	fieldSize = this->getDataSize(1+3*i);
  	partNumber = ms2Data->numberOfParticles[i];
    logger->debug() << "allocating " << partNumber << " floats " << std::endl;
  	logger->debug() << "The field itself reserves place for " << fieldSize << " particles" << std::endl;
  	partStream[i].allocateMemory(partNumber * sizeof(float) * (3 + (colouringVal > 1)) + 1 * sizeof(int));
  	logger->debug() << "inserting number of particles in component "<< i << ": " << partNumber << std::endl;
  	partStream[i] << partNumber;
  	for (int j = 0; j < partNumber; j++)
  	{
  		for (int k = 0; k < 3; k++)
  		{
  			partStream[i] << (float) ((double*) this->getData(1+3*i))[fieldSize*k+j];
  		}
  	}
  }

  // colouringVal = 2: force,
  //              = 3: velocity
  if (colouringVal > 1)
  {
  	int vIndex = colouringVal - 1;
  	logger->debug() << "colouringVal was greater than 1 (" << colouringVal <<") so set some additional data" << std::endl;
  	float tempVal;
  	double* tempPtr;
  	for (int i = 0; i < numberOfComponents; i++)
  	{
  		fieldSize = this->getDataSize(1+3*i);
  		partNumber = ms2Data->numberOfParticles[i];
  		for (int j = 0; j < partNumber; j++)
  		{
  			tempPtr = (double*) this->getData(1 + 3*i + vIndex);

  			tempVal = tempPtr[j]*tempPtr[j] + tempPtr[fieldSize+j]*tempPtr[fieldSize+j]
  			        + tempPtr[2*fieldSize+j]*tempPtr[2*fieldSize+j];
  			tempVal = sqrt(tempVal);
  			if (tempVal > maxVal)
  			{
  				maxVal = tempVal;
  			}
  			if (tempVal < minVal)
  			{
  				minVal = tempVal;
  			}
  			partStream[i] << (float) tempVal;
  		}
  	}
  }
  else
  {
  	minVal = 0;
  	maxVal = (float) numberOfComponents - 1;
  }
  processStream << minVal << maxVal;

  for (int i = 0; i < numberOfComponents; i++)
  {
  	logger->debug() << "inserting partStream[" << i << "] of size " << partStream[i].getStreamSize() << std::endl;
  	processStream << partStream[i];
  	logger->debug() << "size of processStream is now " << processStream.getStreamSize() << std::endl;
  }

  this->getIntraCommunicator()->gatherOnRoot (&processStream, &wholeStream);
  if (this->getIntraCommunicator()->amIRoot ())
  {
  	logger->debug() << "snapshotCommand: assembled the data to be sent in a stream: " << wholeStream->getStreamSize()
    	  						<< std::endl;
  	this->getCommunicator()->sendStream (this->getConnection()->getDataConnection(), *wholeStream);
  }

  logger->debug() << "snapshotCommand: sent the stream" << std::endl;

  if (stepInterval > 0)
  {
    return REPETITION_REQUESTED;
  }
  else
  {
    return EXECUTED;
  }


}

void MegaMolSnapshotCommand::setParameters (std::list <std::string> params)
{
  int listSize = params.size ();
  //std::cout << "set parameters has " << params.size () << " parameters" << std::endl;
  if (listSize > 0)
  {
    stepInterval = atoi (params.front().c_str ());
    params.pop_front();
  }
  if (listSize > 1)
  {
  	colouringVal = atoi (params.front().c_str());
  	params.pop_front();
  }
  //std::cout << "set values: " << stepInterval << " " << colouringVal << std::endl;

}

SteereoCommand* MegaMolSnapshotCommand::generateNewInstance ()
{
	//std::cout << "generate new MegaMolSnapshotCommand " << std::endl;
  SteereoCommand* comm = new MegaMolSnapshotCommand ();
  //startStep = sim -> getSimStep ();
  startStep = 0;
  //((MegaMolSnapshotCommand*) comm) -> setSimData (sim);
  return comm;
}

bool MegaMolSnapshotCommand::condition ()
{
  //int actStep = sd -> getSteps ();
	return true;
	/*if (lastInfo != NULL)
	{
		if (stepInterval > 0)
		{
			unsigned long actStep = lastInfo -> stepNumber;
			return (((actStep - startStep) % stepInterval) == 0);
		}
		else
		{
			return true;
		}
	}*/
}
COMMAND_AUTOREG(MegaMolSnapshotCommand)

#endif
