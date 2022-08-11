
from connect import *
def main(low, high, change, iteration):
    plan = get_current('Plan')
    #suppose only one optimization in the plan
    beamset = get_current('BeamSet')
    #record prescription
    dosetype = beamset.Prescription.PrimaryDosePrescription.PrescriptionType
    dosevalue = beamset.Prescription.PrimaryDosePrescription.DoseValue
    dosevolume = beamset.Prescription.PrimaryDosePrescription.DoseVolume
    roiname = beamset.Prescription.PrimaryDosePrescription.OnStructure.Name
    objectives = plan.PlanOptimizations[0].Objective.ConstituentFunctions  
    count=0

    while count<=iteration:
        flag = 0
        for objective in objectives:
            if objective.ForRegionOfInterest.Type in ['Gtv', 'Ctv', 'Ptv']:
                continue  # if roi belong to targets, skip
            if objective.ForRegionOfInterest.Name in ['PTV']:
                continue  # if roi in ring, skip

            else:
                 if hasattr(objective.DoseFunctionParameters, 'FunctionType'):
                    if objective.DoseFunctionParameters.FunctionType in ['MaxDvh']: 
 
                        if  objective.DoseFunctionParameters.PercentVolume<2:
                            continue
                        else:
                            try:
                                value = objective.FunctionValue.FunctionValue
                            except:
                                continue
                            if value < low:
                                objective.DoseFunctionParameters.PercentVolume -= 1
                                flag = 1
                            elif value>high:
                                objective.DoseFunctionParameters.PercentVolume += 1
                                flag = 1
                    else:
                        try:
                            value = objective.FunctionValue.FunctionValue
                        except:
                            continue
             
                        if value < low:
                            objective.DoseFunctionParameters.DoseLevel -= objective.DoseFunctionParameters.DoseLevel*0.04
                            flag = 1
                        elif value>high:
                            objective.DoseFunctionParameters.DoseLevel += objective.DoseFunctionParameters.DoseLevel*0.04
                            flag = 1
                    count=+1
                 else:
                     continue

        if flag == 0:
            break
        beamset.NormalizeToPrescription(RoiName =roiname, DoseValue = dosevalue, DoseVolume = dosevolume,
         PrescriptionType = dosetype, EvaluateAfterScaling = 1)  #1 present true

#-----------------   ----------------------------------------------------------#
if __name__ == '__main__':
    main(low = 0.0002, high = 0.0012, change = 300, iteration=10)


