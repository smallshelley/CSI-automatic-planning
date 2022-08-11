# This script is going to create an IMRT plan for CSI
# Before running this cript, a local reference point "locref", PTV_brain and PTV_spinal should have been defined.


import wpf
import math
from System.Windows import Application, Window
from connect import *

#--------All structures contain "ctv", "gtv", or "ptv" were recognized as target.  

def changeTargetType(roi,roigeometries):
    for organ in roigeometries:
        organName=organ.OfRoi.Name.lower()
        if organ.PrimaryShape!=None and ('gtv' in organName or 'ctv' in organName or 'ptv' in organName):
           roi[organ.OfRoi.Name].Type = "Ptv"
           roi[organ.OfRoi.Name].OrganData.OrganType = "Target"


#--------Define isocenters, the number of iso is depended on the PTV_spinal length L. If L < 38, two isos were needed, else three.

def createIso(patientmodel,structure_set,examination,Z_size):
    Glocref=structure_set.PoiGeometries['locref']
    Giso=structure_set.RoiGeometries['PTV_spinal'].GetCenterOfRoi()

    #--------The values of couch shifting was integer.

    xs=round(Glocref.Point.x-Giso.x)
    ys=round(Glocref.Point.y-Giso.y)
    zs=round(Glocref.Point.z-Giso.z)

    x0=Glocref.Point.x-xs
    y0=Glocref.Point.y-ys
    z0=Glocref.Point.z-zs

    Giso_brain=structure_set.RoiGeometries['PTV_brain'].GetCenterOfRoi()
    zsb=round(Glocref.Point.z-Giso_brain.z)
    patientmodel.CreatePoi(Examination=examination, Point={ 'x': Glocref.Point.x-xs, 'y': Glocref.Point.y-ys, 'z': Glocref.Point.z-zsb }, Volume=0, Name="iso_brain", Color="Yellow", Type="Isocenter")

    F = 38  #-----the maximum available field size

    if Z_size <= F:
        patientmodel.CreatePoi(Examination=examination, Point={ 'x': x0, 'y': y0, 'z': z0 }, Volume=0, Name="iso_up", Color="Yellow", Type="Isocenter")


    elif F <= Z_size <= 2*F-4:
        patientmodel.CreatePoi(Examination=examination, Point={ 'x': x0, 'y': y0, 'z': z0+int(Z_size/2-F) }, Volume=0, Name="iso_up", Color="Yellow", Type="Isocenter")
        patientmodel.CreatePoi(Examination=examination, Point={ 'x': x0, 'y': y0, 'z': z0-int(Z_size/2-F)-1 }, Volume=0, Name="iso_dn", Color="Yellow", Type="Isocenter")

    else:
        print("The length of PTV is too long")
        exit()

    patientmodel.PointsOfInterest['locref'].Type = "LocalizationPoint"

#--------Create External.  

def createExternal(case):
    retval_0 = case.PatientModel.CreateRoi(Name="External", Color="Green", Type="External", TissueName="", RoiMaterial=None)
    retval_0.CreateExternalGeometry(Examination=examination, ThresholdLevel=-400)


#--------Create Auxiliary Ring.  
          
def CreateRing(case,examination, roi, roiA, roiB, roiname, color, margin1, margin2):
    retval_3 = case.PatientModel.CreateRoi(Name=roiname, Color=color, Type="Control", TissueName=None, RoiMaterial=None)
    retval_3.CreateAlgebraGeometry(Examination=examination, Algorithm="Auto", ExpressionA={ 'Operation': "Union", 'SourceRoiNames': [roiA], 'MarginSettings': 
{ 'Type': "Expand", 'Superior': 0, 'Inferior': 0, 'Anterior': margin1, 'Posterior': margin1, 'Right': margin1, 'Left': margin1 } }, 
ExpressionB={ 'Operation': "Union", 'SourceRoiNames': [roiB], 'MarginSettings': { 'Type': "Expand", 'Superior': margin2, 'Inferior': margin2, 'Anterior': margin2, 'Posterior': margin2, 'Right': margin2, 'Left': margin2 } }, 
ResultOperation="Subtraction", ResultMarginSettings={ 'Type': "Contract", 'Superior': 0, 'Inferior': 0, 'Anterior': 0, 'Posterior': 0, 'Right': 0, 'Left': 0 })


#--------Union PTV_brain and PTV_spinal.  

def OrganAUnionOrganB(case,examination, roi, roiA, roiB, roiname):
    roi[roiname].CreateAlgebraGeometry(Examination=examination, Algorithm="Auto", ExpressionA={ 'Operation': "Union", 'SourceRoiNames': [roiA], 'MarginSettings': 
{ 'Type': "Expand", 'Superior': 0, 'Inferior': 0, 'Anterior': 0, 'Posterior': 0, 'Right': 0, 'Left': 0 } }, 
ExpressionB={ 'Operation': "Union", 'SourceRoiNames': [roiB], 'MarginSettings': { 'Type': "Expand", 'Superior': 0, 'Inferior': 0, 'Anterior': 0, 'Posterior': 0, 'Right': 0, 'Left': 0 } }, 
ResultOperation="Union", ResultMarginSettings={ 'Type': "Expand", 'Superior': 0, 'Inferior': 0, 'Anterior': 0, 'Posterior': 0, 'Right': 0, 'Left': 0 })


#--------CreateOrganNotExist
def CreateOrganNotExist(case,roigeometries,roiname,roi_all_new,roi_all):
 
	if roiname.lower() not in roi_all_new: 
	    try:
                patientmodel.CreateRoi(Name=roiname, Color='Blue', Type="Organ")
	    except:
                print("no such organ !")   


#--------Add plan

def AddPlan(case,examination,poi_all,machine,dose,fraction):
    with CompositeAction('Add Treatment plan'):

         retval_3 = case.AddNewPlan(PlanName="CSI_IMRT", PlannedBy="", Comment="", ExaminationName=examination.Name, AllowDuplicateNames=False)
         retval_3.SetDefaultDoseGrid(VoxelSize={ 'x': 0.3, 'y': 0.3, 'z': 0.3 })
         
         retval_4 = retval_3.AddNewBeamSet(Name="CSI", ExaminationName=examination.Name, MachineName=machine, Modality="Photons", TreatmentTechnique="SMLC", PatientPosition='HeadFirstSupine', 
                 NumberOfFractions=fraction, CreateSetupBeams=False, UseLocalizationPointAsSetupIsocenter=False, Comment="", RbeModelReference=None, EnableDynamicTrackingForVero=False, NewDoseSpecificationPointNames=[],
                 NewDoseSpecificationPoints=[], RespiratoryMotionCompensationTechnique="Disabled", RespiratorySignalSource="Disabled")

         retval_4.AddDosePrescriptionToRoi(RoiName="PTV_all", DoseVolume=95, PrescriptionType="DoseAtVolume", DoseValue=dose, RelativePrescriptionLevel=1, AutoScaleDose=False)

         plan_info = case.QueryPlanInfo(Filter={'Name':"CSI_IMRT"})
         print(list(plan_info))
         plan = case.LoadPlan(PlanInfo=plan_info[0])   

         beam_info = plan.QueryBeamSetInfo( Filter = {'Name': 'CSI'} )
         beam_set = plan.LoadBeamSet( BeamSetInfo = beam_info[0] )

         brain_angles = [60,100,180,260,300]
         spinal_angles = [140,180,220]

         if 'iso_brain' in poi_all:
             iso1 = structure_set.PoiGeometries["iso_brain"]
             
             for i, angle in enumerate(brain_angles):
                 
                 retval_0 = beam_set.CreatePhotonBeam(BeamQualityId=r"6", IsocenterData={ 'Position': { 'x': iso1.Point.x, 'y': iso1.Point.y, 'z': iso1.Point.z }, 'NameOfIsocenterToRef': r"CSI1", 'Name': r"CSI1", 'Color': "98, 184, 234" }, 
             Name=r"Brain"+str(angle), Description=r"", GantryAngle=angle, CouchAngle=0, CouchPitchAngle=0, CouchRollAngle=0, CollimatorAngle=0)

                
         if 'iso_up' in poi_all:
             iso2 = structure_set.PoiGeometries["iso_up"]
             for i, angle in enumerate(spinal_angles):
                 retval_1 = beam_set.CreatePhotonBeam(BeamQualityId=r"6", IsocenterData={ 'Position': { 'x': iso2.Point.x, 'y': iso2.Point.y, 'z': iso2.Point.z }, 'NameOfIsocenterToRef': r"CSI2", 'Name': r"CSI2", 'Color': "98, 184, 234" }, 
             Name=r"Spinal1"+str(angle), Description=r"", GantryAngle=angle, CouchAngle=0, CouchPitchAngle=0, CouchRollAngle=0, CollimatorAngle=0) 


         if 'iso_dn' in poi_all:
             iso3 = structure_set.PoiGeometries["iso_dn"]
             for i, angle in enumerate(spinal_angles):
                 retval_2 = beam_set.CreatePhotonBeam(BeamQualityId=r"6", IsocenterData={ 'Position': { 'x': iso3.Point.x, 'y': iso3.Point.y, 'z': iso3.Point.z }, 'NameOfIsocenterToRef': r"CSI3", 'Name': r"CSI3", 'Color': "98, 184, 234" }, 
             Name=r"Spinal2"+str(angle), Description=r"", GantryAngle=angle, CouchAngle=0, CouchPitchAngle=0, CouchRollAngle=0, CollimatorAngle=0) 


         opt = plan.PlanOptimizations[0]
         opt_param = opt.OptimizationParameters
         opt_param.DoseCalculation.IterationsInPreparationsPhase = 30
         opt_param.Algorithm.MaxNumberOfIterations = 100
         opt_param.DoseCalculation.ComputeFinalDose = True
         opt_param.TreatmentSetupSettings[0].SegmentConversion.MinEquivalentSquare = 4
         opt_param.TreatmentSetupSettings[0].SegmentConversion.MinNumberOfOpenLeafPairs = 2
         opt_param.TreatmentSetupSettings[0].SegmentConversion.MinSegmentArea = 4
         opt_param.TreatmentSetupSettings[0].SegmentConversion.MinSegmentMUPerFraction = 4
         opt_param.TreatmentSetupSettings[0].SegmentConversion.MaxNumberOfSegments = 100

         opt_param.SaveRobustnessParameters(PositionUncertaintyAnterior=0, PositionUncertaintyPosterior=0, PositionUncertaintySuperior=0.5, PositionUncertaintyInferior=0.5, PositionUncertaintyLeft=0, PositionUncertaintyRight=0, 
                   DensityUncertainty=0, PositionUncertaintySetting="IndependentIsocenters", IndependentLeftRight=False, IndependentAnteriorPosterior=False, IndependentSuperiorInferior=True, ComputeExactScenarioDoses=False, 
                   NamesOfNonPlanningExaminations=[])

	 for j in range(5):
             opt_param.TreatmentSetupSettings[0].BeamSettings[j].EditBeamOptimizationSettings(OptimizationTypes=["SegmentOpt", "SegmentMU"], SelectCollimatorAngle=False, AllowBeamSplit=False, JawMotion=r"Use limits as max", 
                   LeftJaw=-15, RightJaw=15, TopJaw=-15, BottomJaw=15)

 # CompositeAction ends 


#--------AddObjectivesTarget

def AddObjectivesTargets(dose):

    MinMaxDoseFuc(opt, 'MinDose', dose, 300, "PTV_all", True)
    DvhFuc(opt, 'MinDvh', dose*1.01, 96, 300, "PTV_all", True)
    MinMaxDoseFuc(opt, 'MaxDose', dose*1.04, 300, "PTV_all", True)


#--------AddObjectivesOAR

def AddObjectivesOARs(structure_set, dose):

    organlist = ['Lung_L','Lung_R','Lung_All','Heart','Lens_L','Lens_R','Kidney_L', 'Kidney_R', 'Kidney_All']
 
    for organ in organlist:
        CreateOrganNotExist(case,roigeometries,organ,roi_all_new,roi_all)

    if roigeometries['Lung_All'].PrimaryShape==None :
        OrganAUnionOrganB(case,examination, roi, 'Lung_L', 'Lung_R', 'Lung_All')
        
    if roigeometries['Kidney_All'].PrimaryShape==None :
        OrganAUnionOrganB(case,examination, roi, 'Kidney_L', 'Kidney_R', 'Kidney_All')

    #--------Get lungs volume
    VolLung = structure_set.RoiGeometries['Lung_All'].GetRoiVolume()

    #--------OARs prediction model was created before

    MaxEudFuc(opt, 'MaxEUD', (0.1*VolLung+446)*dose/2340, 1, 1, 'RingBrain', False)

    MaxEudFuc(opt, 'MaxEUD', (0.124*VolLung+1376)*dose/2340, 1, 1, 'Ring1Spinal', False)

    MaxEudFuc(opt, 'MaxEUD', (0.1*VolLung+890)*dose/2340, 1, 1, 'Ring2Spinal', False)

    MaxEudFuc(opt, 'MaxEUD', (-0.077*VolLung+490)*dose/2340, 1, 1, 'Lung_All', False)

    MaxEudFuc(opt, 'MaxEUD', (-0.053*VolLung+601)*dose/2340, 1, 1, 'Heart', False)

    MaxEudFuc(opt, 'MaxEUD', (-0.070*VolLung+374)*dose/2340, 1, 0.5, 'Kidney_L', False)

    MaxEudFuc(opt, 'MaxEUD', 400, 1, 1, 'Lens_L', False)
    
    MaxEudFuc(opt, 'MaxEUD', 400, 1, 1, 'Lens_R', False)
    
    MinMaxDoseFuc(opt, 'MaxDose', 800, 50, 'Lens_L', False)
    
    MinMaxDoseFuc(opt, 'MaxDose', 800, 50, 'Lens_R', False)

    DoseFallOffFuc(opt, 'DoseFallOff', dose, dose*0.5, 10, 2.5, 'External')


#--------Define Objective Functions


def MinMaxDoseFuc(opt,Function, dose, weight, roiname, isrobust):
    with CompositeAction('Add Optimization Function'):
        retval = opt.AddOptimizationFunction(FunctionType = Function, RoiName = roiname,IsRobust=isrobust)
        retval.DoseFunctionParameters.DoseLevel = dose
        retval.DoseFunctionParameters.Weight = weight


def DoseFallOffFuc(opt, Function, highdose, lowdose, weight, distance, roiname):
    with CompositeAction('Add Optimization Function'):
        retval_f = opt.AddOptimizationFunction(FunctionType = Function, RoiName = roiname)
        retval_f.DoseFunctionParameters.HighDoseLevel =highdose
        retval_f.DoseFunctionParameters.LowDoseLevel =lowdose
        retval_f.DoseFunctionParameters.LowDoseDistance =distance
        retval_f.DoseFunctionParameters.Weight = weight
        
def DvhFuc(opt, Function, dose, percentvolume, weight, roiname, isrobust):
    with CompositeAction('Add Optimization Function'):
        retval_d = opt.AddOptimizationFunction(FunctionType = Function, RoiName = roiname,IsRobust=isrobust)
        retval_d.DoseFunctionParameters.PercentVolume =percentvolume
        retval_d.DoseFunctionParameters.DoseLevel = dose
        retval_d.DoseFunctionParameters.Weight = weight
        
def MaxEudFuc(opt, Function, dose, eudA, weight, roiname, isrobust):
    with CompositeAction('Add Optimization Function'):
        retval_e = opt.AddOptimizationFunction(FunctionType = Function, RoiName = roiname,IsRobust=isrobust)
        retval_e.DoseFunctionParameters.EudParameterA = eudA 
        retval_e.DoseFunctionParameters.DoseLevel = dose
        retval_e.DoseFunctionParameters.Weight = weight  



if __name__ == '__main__':

    patient = get_current("Patient")
    case = get_current("Case")
    examination = get_current("Examination")
    structure_set = case.PatientModel.StructureSets[examination.Name]
    patientmodel = case.PatientModel
    roi = case.PatientModel.RegionsOfInterest
    roigeometries = structure_set.RoiGeometries
    poi = patientmodel.PointsOfInterest
    
    PTV_Bounding = roigeometries['PTV_spinal'].GetBoundingBox() 

    X_size = PTV_Bounding[1].x-PTV_Bounding[0].x
    Y_size = PTV_Bounding[1].y-PTV_Bounding[0].y
    Z_size = PTV_Bounding[1].z-PTV_Bounding[0].z 

   
    createIso(patientmodel,structure_set,examination,Z_size)

    CreateRing(case,examination, roi, "PTV_brain", "PTV_brain", "RingBrain", "Pink", 1.5,  0.5)
    CreateRing(case,examination, roi, "PTV_spinal", "PTV_spinal", "Ring1Spinal", "White", 1.5, 0.5)
    CreateRing(case,examination, roi, "PTV_spinal", "PTV_spinal", "Ring2Spinal", "Yellow", 2.5, 1.5)

    try:
        createExternal(case)
    except:
        pass


    roi_all = [r.Name for r in roi]
    poi_all = [p.Name for p in poi]
    roi_all_new=[]
    for name in roi_all:
        roi_all_new.append(name.lower())

    CreateOrganNotExist(case,roigeometries,'PTV_all',roi_all_new,roi_all)
    OrganAUnionOrganB(case,examination, roi, 'PTV_brain','PTV_spinal','PTV_all')

    changeTargetType(roi,roigeometries)


        
    machine = "3862"  
    PreDose = 2340    
    fraction = 13     

    AddPlan(case, examination, poi_all, machine, PreDose, fraction)
    plan_info=case.QueryPlanInfo(Filter={'Name':"CSI_IMRT"})
    plan = case.LoadPlan(PlanInfo=plan_info[0]) 
    opt = plan.PlanOptimizations[0] 

    AddObjectivesTargets(PreDose)
    AddObjectivesOARs(structure_set, PreDose)

    patient.Save()

    patient.Cases['Case 1'].TreatmentPlans['CSI_IMRT'].PlanOptimizations[0].RunOptimization()
    patient.Cases['Case 1'].TreatmentPlans['CSI_IMRT'].PlanOptimizations[0].RunOptimization()

    for k in range(2):
        execfile(r'I:\wxt\CSI Automatic Planning\SelfAdjust.py')
        patient.Cases['Case 1'].TreatmentPlans['CSI_IMRT'].PlanOptimizations[0].RunOptimization()

    patient.Save()
   
    MinMaxDoseFuc(opt, 'MaxDose', PreDose*1.075, 99999, 'External', False)
    patient.Cases['Case 1'].TreatmentPlans['CSI_IMRT'].PlanOptimizations[0].RunOptimization()

   






