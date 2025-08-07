/*---------------------------------------------------------------------------*\
     ██╗████████╗██╗  ██╗ █████╗  ██████╗ █████╗       ███████╗██╗   ██╗
     ██║╚══██╔══╝██║  ██║██╔══██╗██╔════╝██╔══██╗      ██╔════╝██║   ██║
     ██║   ██║   ███████║███████║██║     ███████║█████╗█████╗  ██║   ██║
     ██║   ██║   ██╔══██║██╔══██║██║     ██╔══██║╚════╝██╔══╝  ╚██╗ ██╔╝
     ██║   ██║   ██║  ██║██║  ██║╚██████╗██║  ██║      ██║      ╚████╔╝
     ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝      ╚═╝       ╚═══╝

 * In real Time Highly Advanced Computational Applications for Finite Volumes
 * Copyright (C) 2017 by the ITHACA-FV authors
-------------------------------------------------------------------------------
License
    This file is part of ITHACA-FV
    ITHACA-FV is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    ITHACA-FV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.
    You should have received a copy of the GNU Lesser General Public License
    along with ITHACA-FV. If not, see <http://www.gnu.org/licenses/>.
Description
    Example of a state and boundary condition reconstruction in 3D heat
    transfer problem using EnKF
SourceFiles
    06enKFwDF_3dIHTP.C
\*---------------------------------------------------------------------------*/

#include <iostream>
#include <fstream>
#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "IOmanip.H"
#include "Time.H"
#include "sequentialIHTP.H"
#include "ITHACAutilities.H"
#include <Eigen/Dense>
#define _USE_MATH_DEFINES
#include <cmath>
#include "Foam2Eigen.H"

#include "MUQ/Modeling/Distributions/Gaussian.h"

#include "muq2ithaca.H"
#include "Fang2017filter_wDF.H"

#include "06enKFwDF_3dIHTP.H"

using namespace SPLINTER;

class TutorialUQ5 : public ITHACAmuq::Fang2017filter_wDF
{
    public:
        explicit TutorialUQ5(int argc, char* argv[], int _Nsamples)
            :
            ITHACAmuq::Fang2017filter_wDF(_Nsamples),
            HTproblem(argc, argv)
        {
            setTime(HTproblem.startTime, HTproblem.deltaTime, HTproblem.endTime);
            setObservationSize(HTproblem.getObservationSize());
            setStateSize(HTproblem.getStateSize());
            setObservationTime(HTproblem.observationStartTimestep, HTproblem.observationDeltaTimesteps);
            HTproblem.setProbe(1, Foam::vector(0.91, 0.02 , 0.55)); // Kabir: Temperature probe.
        }
        inverseHeatTransfer_3D HTproblem;

        //--------------------------------------------------------------------------
        /// Project the state and adds the model error

        void stateProjection()
        {
            Info << "\nState projection start" << endl;
            for (int sampI = 0; sampI < getNumberOfSamples(); sampI++)
            {
                Info << "Sample " << sampI + 1 << ", time = " << getTime() << endl;;
                Eigen::VectorXd newState = HTproblem.projectState( stateEns.getSample(sampI), parameterEns.getSample(sampI),getTime(), getTimeStep(), getTime() + HTproblem.deltaTime, modelErrorDensity);
                stateEns.assignSample(sampI, newState);

                // ################### Kabir: Exporting the ensemble of heat flux weights at each time step in the forecast step, ITHACAoutput/projection/HFW0-HFW99.  0:1:Ntimes
                std::string Weight="heatFlux_weights";
                std::string saveFileName= Weight + std::to_string(sampI);

                std::string folderName = "HFW" + std::to_string(getTimeStep());
                std::string folderPath = "ITHACAoutput/projection/" + folderName;
            
                Eigen::Matrix<double, -1, 1> matrixKabir = parameterEns.getSample(sampI);
                ITHACAstream::exportMatrix(matrixKabir, saveFileName, "eigen", folderPath);
                // ################### Kabir: Exporting the ensemble of heat flux weights at each time step in the forecast step, ITHACAoutput/projection/HFW0-HFW99.   0:1:Ntimes

            }
            Info << "\nState projection end" << endl;
        };

        //--------------------------------------------------------------------------
        /// Observe the state ensamble
        void observeState()
        {
            // Create a matrix 'observedState' with rows based on the size of observations and columns based on the number of samples.
            Eigen::MatrixXd observedState(HTproblem.observe(stateEns.getSample(0)).size(), getNumberOfSamples());                                         // [100 ,300]

            for (int sampI = 0; sampI < getNumberOfSamples(); sampI++)
            {
                // Observe the state for the current sample and store it in the matrix. Add measurement noise to the observed state to make it realistic.
                observedState.col(sampI) = HTproblem.observe(stateEns.getSample(sampI)) + measNoiseDensity->Sample();                                      // [,300]
            }
           // Assign the matrix 'observedState' to an object or variable 'observationEns' for further processing or analysis.
            observationEns.assignSamples(observedState);                                                                                                   // [,300]
        };

        //--------------------------------------------------------------------------
        /// Set parameter(weight) prior mean by projecting initial gTrue(heat flux) on the radial basis functions
              // Kabir: This function is returning parameter(weight) prior mean (in zero time). we need to solve a linear system of equations to get the weights. More information, comment below
              // Kabir: Taking gTrue(the true heat flux which is a function of x in the time zero) and projecting it on the radial basis function. 
                  // Kabir gTrue at the beginning is in the full dimensional space (we have a certain value in each face of the finite volume mesh which is a vector of dimension Nh)
        Eigen::VectorXd setParameterPriorMean()    // q'' = Q'*w --Q--> Qq'' = (Q*Q')*w --> w = (Q*Q')^-1 * Q * q''       We project the true heat flux on the RBFs, It means we start with the correct.
        {
            HTproblem.set_gTrue();
            Eigen::MatrixXd Temp(HTproblem.heatFluxSpaceBasis[0].size(), HTproblem.heatFluxSpaceBasis.size());
            forAll(HTproblem.heatFluxSpaceBasis, baseI)
            {
                Temp.col(baseI) = Foam2Eigen::List2EigenMatrix(HTproblem.heatFluxSpaceBasis[baseI]); // Q'       Kabir: temp[400,5] is equal to the transpose of heatFluxSpaceBasis[25,400], RBF
            }

            Eigen::MatrixXd Btemp = Foam2Eigen::List2EigenMatrix(HTproblem.gTrue[0]);                // q''      Kabir: Btemp[400,1] containing the values of the first row of gTrue (the true heat flux at time zero for all faces of the hotside boundary) 
            Eigen::MatrixXd B = Temp.transpose() * Btemp;                                            // Q * q''  Kabir:[5,1] = [5,400]*[400,1]
            cnpy::save(Btemp, "Btemp.npy");
            cnpy::save(B, "B.npy"); 

            cnpy::save(Temp, "Temp.npy");                                                            //          Kabir: temp here is transpose of the RBF,   RBF <= measurements
            Temp = Temp.transpose() * Temp;                                                          // Q*Q'     Kabir: [5,5] =  [5,400]*[400,5]
            cnpy::save(Temp, "Temp2.npy");                                                           //          Kabir: If this square matrix has everywhere equal numbers it is bad.

            std::cout << "Condition Number of Projection " << EigenFunctions::condNumber(Temp) << std::endl; 
            // Kabir: A condition number close to 1 indicates that the matrix is well-conditioned and small changes in the input will result in small changes in the 
                // output. On the other hand, a large condition number (much greater than 1) indicates that the matrix is ill-conditioned, and small changes in the input can lead to significant changes in the output. 
                         // If I get the condition number of projection that is very high, It means that I get crap in terms of weights. That is because I have
                         // RBFS which are not really orthogonal onto the domain. So, when I am solving the problem associated with the projection , we get really bad weights.So 
                        // if I start with a prior which is really bad, there is no way to get sth better with filtering. 
            
            std::ofstream outputFile("condNumber.txt", std::ios::out);
            outputFile << EigenFunctions::condNumber(Temp);
            outputFile.close();

            //exit(0);
            Eigen::VectorXd result = Temp.fullPivLu().solve(B);                                      // w = (Q*Q')^-1 * Q * q''   Kabir X[5,1]= temp^-1[5,5] * B[5,1] 
            cnpy::save(result, "parameterPriorMeanWithoutShifting.npy"); 
            double ShiftingWeightFactor = 0.3;                 // It must be in the range [0, 1] where 0 corresponds to no modification, and 1 corresponds to increasing each element by 100%.
            // Modify each element of the Projected Prior Weight by adding Shifting Weight Factor % of its current value
            // If ShiftingWeightFactor is set to 0, no modification occurs (i.e., each element remains unchanged), and if ShiftingWeightFactor is set to 1, each element would be doubled.

            for (int i = 0; i < result.size(); ++i)
            {
                result(i) += ShiftingWeightFactor * result(i); // adding ShiftingWeightFactor * result(i) to the current value of result(i).
            }

            return result;                                     // Return the modified Prior Weight 

        };

        //--------------------------------------------------------------------------
        /// Post-processing
        void postProcessing(word outputFolder)
        {
            volScalarField T(HTproblem._T());
            PtrList<volScalarField> TtrueList;
            ITHACAstream::read_fields(TtrueList, "Tdirect","./ITHACAoutput/true/");                                         // TtrueList a shared pointer to read the Tdirect at each time steps

            Eigen::VectorXd probe_rec(getTimeVector().size() - 1); // Kabir:  Temperature probe
            Eigen::VectorXd probeState_maxConf(getTimeVector().size() - 1);
            Eigen::VectorXd probeState_minConf(getTimeVector().size() - 1);

            Eigen::MatrixXd gTrue_probe(1,getTimeVector().size() - 1);
            Eigen::MatrixXd gRec_probe(1,getTimeVector().size() - 1);

            Eigen::MatrixXd gRec_probeMaxConf(1,getTimeVector().size() - 1); // Kabir
            Eigen::MatrixXd gRec_probeMinConf(1,getTimeVector().size() - 1); // Kabir
            

            Foam::vector hotSide_probeLocation(0.91, 0.0 , 0.55);   // Kabir : Heat Flux probe

            for (int timeI = 0; timeI < getTimeVector().size() - 1; timeI++)
            {
                Eigen::VectorXd mean = getStateMean().col(timeI);                                                           // getStateMean().col(timeI) returns stateMean.col(timeI);
                probe_rec(timeI)          = HTproblem.fieldValueAtProbe(mean, HTproblem.probePosition);                     // return output = mean(mesh.findCell(probePosition));
                probeState_maxConf(timeI) = HTproblem.fieldValueAtProbe(state_maxConf.col(timeI), HTproblem.probePosition); // return output = state_maxConf.col(timeI)(mesh.findCell(probePosition));
                probeState_minConf(timeI) = HTproblem.fieldValueAtProbe(state_minConf.col(timeI), HTproblem.probePosition);
                
                volScalarField meanField  = Foam2Eigen::Eigen2field(T, mean);
                ITHACAstream::exportSolution(meanField,       std::to_string(getTime(timeI)), outputFolder,"stateMean");    // stateMean         containing BCs and mean states at each time steps
                ITHACAstream::exportSolution(TtrueList[timeI],std::to_string(getTime(timeI)), outputFolder,"trueState");    // trueState=Tdirect containing BCs and true states at each time steps

                volScalarField diff = TtrueList[timeI] - meanField;
                ITHACAstream::exportSolution(diff, std::to_string(getTime(timeI)), outputFolder, "error");                  // error = trueState=Tdirect - stateMean            at each time steps
                
                volScalarField relativeErrorField(meanField);                                                               // relativeErrorField will contain the same data as meanField
                double EPS = 1e-16;

                for (label i = 0; i < relativeErrorField.internalField().size(); i++) // 3200
                {
                    if (std::abs(TtrueList[timeI].ref()[i]) < EPS)  // TtrueList[timeI].ref()[i] accesses the i-th element of TtrueList[timeI] and std::abs() is used to calculate its absolute value.
                    {
                        relativeErrorField.ref()[i] = (std::abs(diff.ref()[i])) / EPS; 
                    }
                    else
                    {
                        relativeErrorField.ref()[i] = (std::abs(diff.ref()[i])) / TtrueList[timeI].ref()[i];
                    }
                }
                // creates an object of type volScalarField named gTrueField, and then it uses the list2Field function to convert list HTproblem.gTrue[timeI] into a field format, 
                  // which is then assigned to gTrueField. So, gTrueField ends up containing the field data derived from HTproblem.gTrue[timeI]. This allows you to work with the data in 
                  // gTrueField as if it were a field within your code.
                volScalarField gTrueField    = HTproblem.list2Field(HTproblem.gTrue[timeI]);
                ITHACAstream::exportSolution(gTrueField,  std::to_string(HTproblem.timeSteps[timeI]), outputFolder, "gTrue");
                gTrue_probe.col(timeI)       = HTproblem.fieldValueAtProbe(gTrueField, hotSide_probeLocation);                 // return output = gTrueField(mesh.findCell(hotSide_probeLocation));

                volScalarField gField        = HTproblem.list2Field(HTproblem.updateHeatFlux( getParameterMean().col(timeI))); // getParameterMean().col(timeI) returns parameterMean.col(timeStepI);
                ITHACAstream::exportSolution(gField, std::to_string(HTproblem.timeSteps[timeI]), outputFolder, "gRec");
                gRec_probe.col(timeI)        = HTproblem.fieldValueAtProbe(gField,        hotSide_probeLocation);


                volScalarField gFieldMaxConf = HTproblem.list2Field(HTproblem.updateHeatFlux( getParameterMaxConf().col(timeI)));     // Kabir
                ITHACAstream::exportSolution(gFieldMaxConf, std::to_string(HTproblem.timeSteps[timeI]), outputFolder, "gRecMaxConf"); // Kabir
                gRec_probeMaxConf.col(timeI) = HTproblem.fieldValueAtProbe(gFieldMaxConf, hotSide_probeLocation);

                volScalarField gFieldMinConf = HTproblem.list2Field(HTproblem.updateHeatFlux( getParameterMinConf().col(timeI)));     // Kabir
                ITHACAstream::exportSolution(gFieldMinConf, std::to_string(HTproblem.timeSteps[timeI]), outputFolder, "gRecMinConf"); // Kabir
                gRec_probeMinConf.col(timeI) = HTproblem.fieldValueAtProbe(gFieldMinConf, hotSide_probeLocation);

                ITHACAstream::exportSolution(relativeErrorField, std::to_string(getTime(timeI)), outputFolder, "relativeErrorField"); // relativeErrorField                    at each time steps
            }
            ITHACAstream::exportMatrix(probe_rec, "probe_rec", "eigen", outputFolder);
            ITHACAstream::exportMatrix(probeState_maxConf, "probeState_maxConf", "eigen", outputFolder);
            ITHACAstream::exportMatrix(probeState_minConf, "probeState_minConf", "eigen", outputFolder);

            ITHACAstream::exportMatrix(gTrue_probe, "gTrue_probe", "eigen", outputFolder);
            ITHACAstream::exportMatrix(gRec_probe, "gRec_probe", "eigen", outputFolder);

            ITHACAstream::exportMatrix(gRec_probeMaxConf, "gRec_probeMaxConf", "eigen", outputFolder); // Kabir
            ITHACAstream::exportMatrix(gRec_probeMinConf, "gRec_probeMinConf", "eigen", outputFolder); // Kabir
        }

};

int main(int argc, char* argv[])
{
    int Nsamples = 300; // Kabir: All the ensemble-based methods, in general, tend to converge as we increase the number of samples. Therefore, the number of samples should be as high as possible (the more sample we have usually the more accurate).
    
    TutorialUQ5 example(argc, argv, Nsamples);

    // Reading parameters from file
    ITHACAparameters* para = ITHACAparameters::getInstance( example.HTproblem._mesh(), example.HTproblem._runTime());
    
    example.HTproblem.a = para->ITHACAdict->lookupOrDefault<scalar>("a", 0);
    example.HTproblem.b = para->ITHACAdict->lookupOrDefault<scalar>("b", 0);
    example.HTproblem.c = para->ITHACAdict->lookupOrDefault<scalar>("c", 0);
    example.HTproblem.maxFrequency = para->ITHACAdict->lookupOrDefault<scalar>("maxFrequency", 0);

    example.HTproblem.HTC = para->ITHACAdict->lookupOrDefault<scalar>("heatTranferCoeff", 0);
    example.HTproblem.thermalCond = para->ITHACAdict->lookupOrDefault<scalar>("thermalConductivity", 0.0);
    example.HTproblem.density = para->ITHACAdict->lookupOrDefault<scalar>("density", 0.0);
    example.HTproblem.specificHeat = para->ITHACAdict->lookupOrDefault<scalar>("specificHeat", 0.0);
    example.HTproblem.initialField = para->ITHACAdict->lookupOrDefault<scalar>("initialField", 0);

    scalar measNoiseCov = para->ITHACAdict->lookupOrDefault<scalar>("measNoiseCov", 0);
    double modelErrorCov = para->ITHACAdict->lookupOrDefault<double>("modelErrorCov", 0);
    scalar stateCov = para->ITHACAdict->lookupOrDefault<scalar>("stateInitialCov", 0);
    scalar parameterCov = para->ITHACAdict->lookupOrDefault<scalar>("parameterPriorCov", 0);            // Kabir: Inside the code, it is calculated by taking a (20%) percent of parameterPriorMean
    
    label NheatFluxPODbasis = para->ITHACAdict->lookupOrDefault<label>("NheatFluxPODbasis", 0);

    label sizeOfParameter = para->ITHACAdict->lookupOrDefault<label>("sizeOfTheParameter", 0);          // label sizeOfParameter = 5; 
    label innerLoops = para->ITHACAdict->lookupOrDefault<label>("EnKF_innerLoop", 1);
    scalar basisShapeParameter = para->ITHACAdict->lookupOrDefault<scalar>("basisofShapeParameter", 0); // scalar basisShapeParameter; 
    word reconstructionFolder = "ITHACAoutput/reconstruction";

    example.HTproblem.setSpaceBasis("rbf", basisShapeParameter, NheatFluxPODbasis);                     /// Define the RBFs used for the parametrization of g
    example.HTproblem.Nbasis = sizeOfParameter; 

    const int stateSize = example.getStateSize();

    example.setParameterSize(sizeOfParameter);                                                    // Kabir: parametersize is initialized as sizeOfParameter = 5
    const int parameterSize = example.getParameterSize();                                         // Kabir: return parametersize = 5

    Eigen::VectorXd stateInitialMean =Eigen::VectorXd::Ones(1) * example.HTproblem.initialField;  // stateInitialMean= initialField = 400
    Eigen::MatrixXd stateInitialCov = Eigen::MatrixXd::Identity(1,1) * stateCov;                  // stateInitialCov = stateCov     = 10

    Eigen::VectorXd parameterPriorMean = example.setParameterPriorMean();                         // Kabir: Please see comments in front of setParameterPriorMean function to understand how it works. The size of the parameterPriorMean is equal to sizeOfParameter = 5
    
    // double scaleFactor = 0.2;   
    double scaleFactor = para->ITHACAdict->lookupOrDefault<double>("scaleOfFactor", 0);           // Kabir: For the covariance of the prior weights, we can take, for example, 20 percent of the weight prior mean(parameterPriorMean) by defining scaleFactor
    double smallNumber = 1e-3;                                                                    // Kabir: We can adjust the value of the small number as needed
    Eigen::VectorXd parameterPriorMean1 = parameterPriorMean;                                     // Kabir: Create a copy of parameterPriorMean
               // Addinng a really small number to all the zero diagonal elements of parameterPriorCov if parameterPriorMean vector has probaboly some zero elements.(parameterPriorMean should not have zero element)
    for (int i = 0; i < parameterPriorMean.rows(); ++i) {
        if (parameterPriorMean1(i) == 0) {
            parameterPriorMean1(i) += smallNumber;
        }
    }
    std::string WeightPriorMean1="prior_weights_Mean1";
                                                                                                   // Compute the diagonal covariance matrix parameterPriorCov
    Eigen::MatrixXd parameterPriorCov = scaleFactor * parameterPriorMean1.cwiseAbs().asDiagonal(); // creates a diagonal matrix parameterPriorCov where each diagonal element is obtained by scaling the corresponding element of parameterPriorMean by the scaleFactor
                                                                                                   // Eigen::VectorXd class does not have a member function named abs, I used cwiseAbs() but got the same error.

    // ################### Kabir: Exporting parameterPriorMean and parameterPriorCov in order to plot the prior PDF of parameters.
    std::string WeightPriorMean="prior_weights_Mean";
    std::string WeightPriorCov="prior_weights_Cov";
    std::string folderNamePr = "PriorMeanCovariance";
    std::string folderPathPr = "ITHACAoutput/projection/" + folderNamePr;
                                                                                                   //Eigen::MatrixXd MatrixPriroCov = example.parameterPriorDensity->GetCovariance();  // Kabir: No need, but Giovanni told me use this line to GetCovariance, source from this website https://mituq.bitbucket.io/source/_site/latest/classmuq_1_1Modeling_1_1Gaussian.html
    ITHACAstream::exportMatrix(parameterPriorMean, WeightPriorMean, "eigen", folderPathPr);
    ITHACAstream::exportMatrix(parameterPriorCov, WeightPriorCov, "eigen", folderPathPr);
    cnpy::save(parameterPriorMean, "parameterPriorMean.npy");                                      // Save parameterPriorMean and parameterPriorCov as a numpy array
    cnpy::save(parameterPriorCov, "parameterPriorCov.npy");
    ITHACAstream::exportMatrix(parameterPriorMean1, WeightPriorMean1, "eigen", folderPathPr);
    // ################### Kabir: Exporting the mean vector and covariance matrix of weights in order to plot the prior PDF of parameters.
    
    Eigen::MatrixXd measurementsMat = example.HTproblem.solveDirect();                             // [25,50] 
    cnpy::save(measurementsMat, "measurementsMat.npy");  
    ITHACAstream::exportMatrix(measurementsMat, "measurementsMat_noNoise", "eigen", reconstructionFolder); 

    // Add noise to measurements, as a variance with respect to what we are reading
    example.setMeasNoise(measNoiseCov * measurementsMat.mean());                                   // Kabir: Covariance measurement noise = measNoiseCov(0.0001) * measurementsMat.mean()(340.5880) =0.0340588, as a variance with respect to what we are reading,
    // ################### Kabir:
    double measNoiseCovDotmeasurementsMatMean = measNoiseCov * measurementsMat.mean();             //Kabir: measurementsMat.mean() = (sum of all elements) / (M * N)
    std::cout << "Name: measNoiseCov * measurementsMat.mean() is equal to " << std::endl;
    std::cout << measNoiseCovDotmeasurementsMatMean << std::endl;                                  // 0.0340588
    // Open a file for output
    std::ofstream outputFile("measNoiseCovTotal.txt", std::ios::out);
    // Write the value of measNoiseCovDotmeasurementsMatMean to the file
    outputFile << measNoiseCovDotmeasurementsMatMean;
    outputFile.close();

    // ################### Kabir:
    for(int i = 0; i < measurementsMat.cols(); i++)                                                // [25,50]
    {
        measurementsMat.col(i) = measurementsMat.col(i) + example.measNoiseDensity->Sample(); 
    }
    ITHACAstream::exportMatrix(measurementsMat, "measurementsMat_noise", "eigen", reconstructionFolder);
    cnpy::save(measurementsMat, "measurementsMatNoise.npy");

    example.setObservations(measurementsMat);

    bool univariateInitStateDensFlag = 1;
    example.setInitialStateDensity(stateInitialMean, stateInitialCov, univariateInitStateDensFlag);
    std::cout << "debug: parameterPriorMean = " << parameterPriorMean << std::endl;
    example.setParameterPriorDensity(parameterPriorMean, parameterPriorCov);

    bool univariateModelErrorDistribution = 1;
    example.setModelError(modelErrorCov, univariateModelErrorDistribution);
    //example.setMeasNoise(0.1);                                                                    // Kabir: No need, because, some lines above, we are using example.setMeasNoise(measNoiseCov * measurementsMat.mean());
    example.run(innerLoops, reconstructionFolder);
    example.postProcessing(reconstructionFolder);
    return 0;
}
