//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B2aDetectorConstruction.cc
/// \brief Implementation of the B2aDetectorConstruction class
 
#include "B2aDetectorConstruction.hh"
#include "B2aDetectorMessenger.hh"
#include "B2TrackerSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4ThreadLocal 
G4GlobalMagFieldMessenger* B2aDetectorConstruction::fMagFieldMessenger = 0;

B2aDetectorConstruction::B2aDetectorConstruction()
:G4VUserDetectorConstruction(), 
 fNbOfChambers(0),
 fLogicTarget(NULL), fLogicChamber(NULL), 
 fTargetMaterial(NULL), fChamberMaterial(NULL), 
 fStepLimit(NULL),
 fCheckOverlaps(true)
{
  fMessenger = new B2aDetectorMessenger(this);

  fNbOfChambers = 5;
  fLogicChamber = new G4LogicalVolume*[fNbOfChambers];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
B2aDetectorConstruction::~B2aDetectorConstruction()
{
  delete [] fLogicChamber; 
  delete fStepLimit;
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4VPhysicalVolume* B2aDetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2aDetectorConstruction::DefineMaterials()
{
  // Material definition 

  G4NistManager* nistManager = G4NistManager::Instance();

  // Air defined using NIST Manager
  nistManager->FindOrBuildMaterial("G4_AIR");
  
  // Lead defined using NIST Manager
  fTargetMaterial  = nistManager->FindOrBuildMaterial("G4_Pb");

  // Xenon gas defined using NIST Manager
  fChamberMaterial = nistManager->FindOrBuildMaterial("G4_Xe");

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B2aDetectorConstruction::DefineVolumes()
{

  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  
  //Material stuff

    G4double a;  // mass of a mole;
    G4double z;  // z=mean number of protons;
    G4double density, fractionMass;
    G4String symbol, name;
    G4int nComponents, nAtoms;
    //G4double temp;

    G4Element* elH  = new G4Element(name = "Hydrogen"   , symbol = "H"  , z = 1.  , a =  1.008*g/mole);
    G4Element* elHe = new G4Element(name = "Helium"     , symbol = "He" , z = 2.  , a =  4.003*g/mole);
    G4Element* elC  = new G4Element(name = "Carbon"     , symbol = "C"  , z = 6.  , a = 12.011*g/mole);
    G4Element* elN  = new G4Element(name = "Nitrogen"   , symbol = "N"  , z = 7.  , a = 14.007*g/mole);
    G4Element* elO  = new G4Element(name = "Oxygen"     , symbol = "O"  , z = 8.  , a = 15.999*g/mole);
    G4Element* elSi = new G4Element(name = "Silicon"    , symbol = "Si" , z = 14. , a = 28.086*g/mole);
    G4Element* elP  = new G4Element(name = "Phosphorus" , symbol = "P"  , z = 15. , a = 30.974*g/mole);
    G4Element* elS  = new G4Element(name = "Sulfur"     , symbol = "S"  , z = 16. , a = 32.065*g/mole);
    G4Element* elCr = new G4Element(name = "Chromium"   , symbol = "Cr" , z = 24. , a = 51.996*g/mole);
    G4Element* elMn = new G4Element(name = "Manganese"  , symbol = "Mn" , z = 25. , a = 54.938*g/mole);
    G4Element* elFe = new G4Element(name = "Iron"       , symbol = "Fe" , z = 26. , a = 55.845*g/mole);
    G4Element* elNi = new G4Element(name = "Nickel"     , symbol = "Ni" , z = 28. , a = 58.693*g/mole);
    G4Element* elMo = new G4Element(name = "Molybdenum" , symbol = "Mo" , z = 42. , a = 95.94 *g/mole);
    G4Element* elAl = new G4Element(name = "Aluminuim"  , symbol = "Al" , z = 13. , a = 26.98 *g/mole);
    G4Element* elMg = new G4Element(name = "Magnesium"  , symbol = "Mg" , z = 12. , a = 24.305*g/mole);
    G4Element* elCu = new G4Element(name = "Copper"     , symbol = "Cu" , z = 29. , a = 63.546*g/mole);
    G4Element* elZn = new G4Element(name = "Zinc"       , symbol = "Zn" , z = 30. , a = 65.38 *g/mole);
    G4Element* elTi = new G4Element(name = "Titanium"   , symbol = "Ti" , z = 22. , a = 47.867*g/mole);
    
    G4Material* liquid_helium = new G4Material(name = "liquid_helium", density= 0.141*g/cm3, nComponents=1); //, kStateLiquid, temp=3.*kelvin);
    liquid_helium -> AddElement(elHe, fractionMass = 100*perCent);
    
    G4Material* liquid_nitrogen = new G4Material(name = "liquid_nitrogen", density= 0.807*g/cm3, nComponents=1); //, kStateLiquid, temp=77.*kelvin);
    liquid_nitrogen -> AddElement(elN, fractionMass = 100*perCent);
    
    //Stainless Steel Type 304
    G4Material* ss_t316 = new G4Material(name = "ss_t316" , density = 8.03*g/cm3 , nComponents = 10);
    ss_t316 -> AddElement(elC  , fractionMass =   0.08*perCent);
    ss_t316 -> AddElement(elMn , fractionMass =   2.00*perCent);
    ss_t316 -> AddElement(elP  , fractionMass =  0.045*perCent);
    ss_t316 -> AddElement(elS  , fractionMass =  0.030*perCent);
    ss_t316 -> AddElement(elSi , fractionMass =   0.75*perCent);
    ss_t316 -> AddElement(elCr , fractionMass =  17.00*perCent);
    ss_t316 -> AddElement(elNi , fractionMass =  12.00*perCent);
    ss_t316 -> AddElement(elMo , fractionMass =   2.50*perCent);
    ss_t316 -> AddElement(elN  , fractionMass =   0.10*perCent);
    ss_t316 -> AddElement(elFe , fractionMass = 65.495*perCent);
    
    //6061 Aluminuim
    G4Material* al_6061 = new G4Material(name = "al_6061" , density = 2.70*g/cm3 , nComponents = 9);
    al_6061 -> AddElement(elAl , fractionMass = 97.3*perCent);
    al_6061 -> AddElement(elMg , fractionMass = 1.0*perCent);
    al_6061 -> AddElement(elSi , fractionMass = 0.6*perCent);
    al_6061 -> AddElement(elFe , fractionMass = 0.3*perCent);
    al_6061 -> AddElement(elCu , fractionMass = 0.3*perCent);
    al_6061 -> AddElement(elCr , fractionMass = 0.2*perCent);
    al_6061 -> AddElement(elZn , fractionMass = 0.1*perCent);
    al_6061 -> AddElement(elTi , fractionMass = 0.1*perCent);
    al_6061 -> AddElement(elMn , fractionMass = 0.1*perCent);

    G4Material* Synthetic_Silica = new G4Material(name = "Synthetic_Silica" , density = 2.65*g/cm3 , nComponents = 2);
    Synthetic_Silica -> AddElement(elSi , nAtoms = 1);
    Synthetic_Silica -> AddElement(elO  , nAtoms = 2);

    G4Material* sapphire = new G4Material(name = "sapphire" , density = 3.98*g/cm3 , nComponents = 2);
    sapphire -> AddElement(elAl , nAtoms = 2);
    sapphire -> AddElement(elO  , nAtoms = 3);

    G4Material* ti_material = new G4Material(name = "ti_material", density = 4.11*g/cm3, nComponents = 1);
    ti_material -> AddElement(elTi, nAtoms = 1);


    G4Material* cu_material = new G4Material(name = "cu_material", density = 8.96*g/cm3, nComponents = 1);
    cu_material -> AddElement(elCu, nAtoms = 1);

  G4Material* liquid_helium_material = G4Material::GetMaterial("liquid_helium");
  G4Material* liquid_nitrogen_material = G4Material::GetMaterial("liquid_nitrogen");
  G4Material* stainless_material = G4Material::GetMaterial("ss_t316");
  G4Material* al_6061_material = G4Material::GetMaterial("al_6061");
  G4Material* silica_material = G4Material::GetMaterial("Synthetic_Silica");
  G4Material* sapphire_material = G4Material::GetMaterial("sapphire");
  
  G4VisAttributes* helium_vis = new G4VisAttributes(G4Colour(0.6,0.0,1.0,0.75));
  G4VisAttributes* nitrogen_vis = new G4VisAttributes(G4Colour(0.0,1.0,1.0,0.75));
  G4VisAttributes* container_vis = new G4VisAttributes(G4Colour(0.5,0.5,0.5,0.75));
  container_vis->SetVisibility(true);

  //===============  Visualization ===============//

  G4VisAttributes* yellowTVA = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0, 0.5));
  G4VisAttributes* redTVA = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0, 0.5));
  G4VisAttributes* greenTVA = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0, 0.5));
  G4VisAttributes* greyTVA = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 0.25));
  G4VisAttributes* blueTVA = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.5));


  // Envelope parameters
  //
  G4double env_sizeXY = 300*cm, env_sizeZ = 300*cm;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_AIR");
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = false;

  
  //     
  // World
  //
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  //G4Material* world_mat = liquid_helium_material;
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name

  logicWorld->SetVisAttributes(greyTVA);
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                     
  //     
  // Envelope
  //  
  G4Box* solidEnv =    
    new G4Box("Envelope",                    //its name
        0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size
      
  G4LogicalVolume* logicEnv =                         
    new G4LogicalVolume(solidEnv,            //its solid
                        env_mat,             //its material
                        "Envelope");         //its name

  logicEnv->SetVisAttributes(greyTVA);
               
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
 

  G4RotationMatrix* yRot = new G4RotationMatrix;  // Rotates X and Z axes only
  yRot->rotateY(90*deg);                     // Rotates 90 degrees

  //Helium Detector
  G4double he_det_rad = 1.875*2.54*cm/2.0;
  G4double he_det_height = 1.743*2.54*cm;
  
  G4ThreeVector he_det_pos = G4ThreeVector(0*cm, 0*cm, 0*cm);
  G4Tubs* he_det_shape =    
     new G4Tubs("he_det", 
         0*cm, 
         he_det_rad, 
         he_det_height/2.0, 
         0.*deg, 
         360.*deg);
  
  G4LogicalVolume* he_det_logic = 
     new G4LogicalVolume(he_det_shape,
                         liquid_helium_material,
                         "he_det");
  
  he_det_logic->SetVisAttributes(helium_vis);
  
  new G4PVPlacement(yRot,                       //no rotation
                    he_det_pos,              //at position
                    he_det_logic,            //its logical volume
                    "he_det",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  
  //Aluminuim Reflector
  G4double al_reflector_OD = 1.98375*2.54*cm;
  G4double al_reflector_ID = 2.0*he_det_rad;
  G4double al_reflector_height = he_det_height;
  
  G4ThreeVector al_reflector_pos = G4ThreeVector(0*cm, 0*cm, 0*cm) + he_det_pos;
  G4Tubs* al_reflector_shape =    
     new G4Tubs("al_reflector", 
         al_reflector_ID/2.0, 
         al_reflector_OD/2.0, 
         al_reflector_height/2.0, 
         0.*deg, 
         360.*deg);
  
  G4LogicalVolume* al_reflector_logic = 
     new G4LogicalVolume(al_reflector_shape,
                         al_6061_material,
                         "al_reflector");
  
  al_reflector_logic->SetVisAttributes(redTVA);
  
  new G4PVPlacement(yRot,                       //no rotation
                    al_reflector_pos,        //at position
                    al_reflector_logic,      //its logical volume
                    "al_reflector",          //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
 
  //SiO2 Windows
  G4double silica_window_OD = al_reflector_OD;
  G4double silica_window_ID = 0*cm;
  G4double silica_window_height = 0.059 *2.54*cm;
  
  G4ThreeVector silica_window_pos_1 =  G4ThreeVector(he_det_height/2.0 + silica_window_height/2.0, 0*cm, 0*cm) + he_det_pos;
  G4ThreeVector silica_window_pos_2 = -G4ThreeVector(he_det_height/2.0 + silica_window_height/2.0, 0*cm, 0*cm) + he_det_pos;
  G4Tubs* silica_window_shape =    
     new G4Tubs("silica_window", 
         silica_window_ID/2.0, 
         silica_window_OD/2.0, 
         silica_window_height/2.0, 
         0.*deg, 
         360.*deg);
  
  G4LogicalVolume* silica_window_logic = 
     new G4LogicalVolume(silica_window_shape,
                         silica_material,
                         "silica_window");
  
  silica_window_logic->SetVisAttributes(redTVA);
  
  new G4PVPlacement(yRot,                       //no rotation
                    silica_window_pos_1,     //at position
                    silica_window_logic,     //its logical volume
                    "silica_window_1",         //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  new G4PVPlacement(yRot,                       //no rotation
                    silica_window_pos_2,     //at position
                    silica_window_logic,     //its logical volume
                    "silica_window_2",         //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  //Sapphire Windows
  G4double sapphire_window_OD = 2.0*2.54*cm;
  G4double sapphire_window_ID = 0*cm;
  G4double sapphire_window_height = 0.01*2.54 * cm;
  
  G4ThreeVector sapphire_window_pos_1 =  G4ThreeVector(0.993*2.54*cm - sapphire_window_height/2.0, 0*cm, 0*cm) + he_det_pos;
  G4ThreeVector sapphire_window_pos_2 = -G4ThreeVector(0.993*2.54*cm - sapphire_window_height/2.0, 0*cm, 0*cm) + he_det_pos;
  G4Tubs* sapphire_window_shape =    
     new G4Tubs("sapphire_window", 
         sapphire_window_ID/2.0, 
         sapphire_window_OD/2.0, 
         sapphire_window_height/2.0, 
         0.*deg, 
         360.*deg);
  
  G4LogicalVolume* sapphire_window_logic = 
     new G4LogicalVolume(sapphire_window_shape,
                         sapphire_material,
                         "sapphire_window");
  
  sapphire_window_logic->SetVisAttributes(greenTVA);
  
  new G4PVPlacement(yRot,                       //no rotation
                    sapphire_window_pos_1,     //at position
                    sapphire_window_logic,     //its logical volume
                    "sapphire_window_1",         //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    false);          //overlaps checking

  new G4PVPlacement(yRot,                       //no rotation
                    sapphire_window_pos_2,     //at position
                    sapphire_window_logic,     //its logical volume
                    "sapphire_window_2",         //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    false);          //overlaps checking

  //Ti flange tube
  G4double ti_flange_tube_OD = 2.050*2.54*cm;
  G4double ti_flange_tube_ID = 1.990*2.55*cm;
  G4double ti_flange_tube_height = 0.572*2.54*cm;
  
  G4ThreeVector ti_flange_tube_pos_1 = G4ThreeVector( 0.425*2.54*cm + ti_flange_tube_height/2.0, 0*cm, 0*cm) + he_det_pos;
  G4ThreeVector ti_flange_tube_pos_2 = G4ThreeVector(-0.425*2.54*cm - ti_flange_tube_height/2.0, 0*cm, 0*cm) + he_det_pos;
  G4Tubs* ti_flange_tube_shape =    
     new G4Tubs("ti_flange_tube", 
         ti_flange_tube_ID/2.0, 
         ti_flange_tube_OD/2.0, 
         ti_flange_tube_height/2.0, 
         0.*deg, 
         360.*deg);
  
  G4LogicalVolume* ti_flange_tube_logic = 
     new G4LogicalVolume(ti_flange_tube_shape,
                         ti_material,
                         "ti_flange_tube");
  
  ti_flange_tube_logic->SetVisAttributes(yellowTVA);
 
  new G4PVPlacement(yRot,                       //no rotation
                    ti_flange_tube_pos_1,    //at position
                    ti_flange_tube_logic,    //its logical volume
                    "ti_flange_tube",        //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    false);          //overlaps checking
  

  new G4PVPlacement(yRot,                       //no rotation
                    ti_flange_tube_pos_2,    //at position
                    ti_flange_tube_logic,    //its logical volume
                    "ti_flange_tube",        //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    false);          //overlaps checking

  //Ti flange inner
  G4double ti_flange_inner_OD = 2.70*2.54*cm;
  G4double ti_flange_inner_ID = ti_flange_tube_ID;
  G4double ti_flange_inner_height = 0.212*2.54*cm;
  
  G4ThreeVector ti_flange_inner_pos_1 = G4ThreeVector( 0.425*2.54*cm + ti_flange_inner_height/2.0, 0*cm, 0*cm) + he_det_pos;
  G4ThreeVector ti_flange_inner_pos_2 = G4ThreeVector(-0.425*2.54*cm - ti_flange_inner_height/2.0, 0*cm, 0*cm) + he_det_pos;
  G4Tubs* ti_flange_inner_shape =    
     new G4Tubs("ti_flange_inner", 
         ti_flange_inner_ID/2.0, 
         ti_flange_inner_OD/2.0, 
         ti_flange_inner_height/2.0, 
         0.*deg, 
         360.*deg);
  
  G4LogicalVolume* ti_flange_inner_logic = 
     new G4LogicalVolume(ti_flange_inner_shape,
                         ti_material,
                         "ti_flange_inner");
  
  ti_flange_inner_logic->SetVisAttributes(yellowTVA);
  
  new G4PVPlacement(yRot,                       //no rotation
                    ti_flange_inner_pos_1,    //at position
                    ti_flange_inner_logic,    //its logical volume
                    "ti_flange_inner",        //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  

  new G4PVPlacement(yRot,                       //no rotation
                    ti_flange_inner_pos_2,    //at position
                    ti_flange_inner_logic,    //its logical volume
                    "ti_flange_inner",        //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking



  //Ti flange outer
  G4double ti_flange_outer_OD = 3.40*2.54*cm;
  G4double ti_flange_outer_ID = ti_flange_tube_OD;
  G4double ti_flange_outer_height = 0.132*2.54*cm;
  
  G4ThreeVector ti_flange_outer_pos_1 = G4ThreeVector( 0.505*2.54*cm + ti_flange_outer_height/2.0, 0*cm, 0*cm) + he_det_pos;
  G4ThreeVector ti_flange_outer_pos_2 = G4ThreeVector(-0.505*2.54*cm - ti_flange_outer_height/2.0, 0*cm, 0*cm) + he_det_pos;
  G4Tubs* ti_flange_outer_shape =    
     new G4Tubs("ti_flange_outer", 
         ti_flange_outer_ID/2.0, 
         ti_flange_outer_OD/2.0, 
         ti_flange_outer_height/2.0, 
         0.*deg, 
         360.*deg);
  
  G4LogicalVolume* ti_flange_outer_logic = 
     new G4LogicalVolume(ti_flange_outer_shape,
                         ti_material,
                         "ti_flange_outer");
  
  ti_flange_outer_logic->SetVisAttributes(yellowTVA);
  
  new G4PVPlacement(yRot,                       //no rotation
                    ti_flange_outer_pos_1,   //at position
                    ti_flange_outer_logic,   //its logical volume
                    "ti_flange_outer",       //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  

  new G4PVPlacement(yRot,                       //no rotation
                    ti_flange_outer_pos_2,   //at position
                    ti_flange_outer_logic,   //its logical volume
                    "ti_flange_outer",       //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  //SS Body
  G4double ss_body_side1 = 3.40*2.54*cm;
  G4double ss_body_side2 = ss_body_side1;
  G4double ss_body_height = 1.01*2.54*cm;

  G4double ss_main_hole_OD = 2.0*2.54*cm;
  G4double ss_main_hole_ID = 0.0*cm;
  G4double ss_main_hole_height = ss_body_height + 1*cm;

  G4double ss_pocket_OD = 2.705*2.54*cm;
  G4double ss_pocket_ID = ss_main_hole_OD - 1*cm;
  G4double ss_pocket_height = 0.105*2.54*cm;
  
  G4ThreeVector ss_pocket_pos_1 = G4ThreeVector(0*cm, 0*cm, ss_body_height/2.0);
  G4ThreeVector ss_pocket_pos_2 = G4ThreeVector(0*cm, 0*cm, -ss_body_height/2.0);
  G4Box* ss_body_shape =    
     new G4Box("ss_body", 
         ss_body_side1/2.0, 
         ss_body_side2/2.0, 
         ss_body_height/2.0);

  G4Tubs* ss_main_hole_shape =    
     new G4Tubs("ss_main_hole", 
         ss_main_hole_ID/2.0, 
         ss_main_hole_OD/2.0, 
         ss_main_hole_height/2.0, 
         0.*deg, 
         360.*deg); 

  G4Tubs* ss_pocket_shape =    
     new G4Tubs("ss_pocket", 
         ss_pocket_ID/2.0, 
         ss_pocket_OD/2.0, 
         ss_pocket_height, 
         0.*deg, 
         360.*deg);

  G4SubtractionSolid* subtracted1 = 
     new G4SubtractionSolid("subtracted1", ss_body_shape, ss_main_hole_shape);
  G4SubtractionSolid* subtracted2 = 
     new G4SubtractionSolid("subtracted2", subtracted1, ss_pocket_shape, 0, ss_pocket_pos_1);
  G4SubtractionSolid* ss_block_shape = 
     new G4SubtractionSolid("ss_block_shape", subtracted2, ss_pocket_shape, 0, ss_pocket_pos_2);

  G4LogicalVolume* ss_block_logic = 
     new G4LogicalVolume(ss_block_shape,
                         ss_t316,
                         "ss_block");
  
  ss_block_logic->SetVisAttributes(greenTVA);
  
  new G4PVPlacement(yRot,                       //no rotation
                    he_det_pos,              //at position
                    ss_block_logic,          //its logical volume
                    "ss_block",              //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    false);          //overlaps checking
  
  //
  //Dewar construction
  //

  //IVC can
  G4double ivc_can_OD = 13.180*2.54*cm;
  G4double ivc_can_thick = 0.090*2.54*cm;
  G4double ivc_can_height = 45.778*2.54*cm;

  G4double dewar_top_to_MC = 39.077*2.54*cm;
  G4double he_det_center_to_MC = 4.0*2.54*cm;
  G4double ivc_can_center_to_dewar_top = 47.464*2.54*cm;
  
  G4ThreeVector ivc_can_pos = G4ThreeVector(0*cm, 0*cm, -ivc_can_center_to_dewar_top + dewar_top_to_MC + he_det_center_to_MC);
  G4Tubs* ivc_can_shape =    
     new G4Tubs("ivc_can", 
         ivc_can_OD/2.0 - ivc_can_thick, 
         ivc_can_OD/2.0, 
         ivc_can_height/2.0, 
         0.*deg, 
         360.*deg);
  
  G4LogicalVolume* ivc_can_logic = 
     new G4LogicalVolume(ivc_can_shape,
                         ss_t316,
                         "ivc_can");
  
  ivc_can_logic->SetVisAttributes(redTVA);
  
  new G4PVPlacement(0,                       //no rotation
                    ivc_can_pos,             //at position
                    ivc_can_logic,           //its logical volume
                    "ivc_can",      //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking


  //dewar inner vacuum can
  G4double dewar_inner_vac_OD = 15.942*2.54*cm;
  G4double dewar_inner_vac_thick = 0.028*2.54*cm;
  G4double dewar_inner_vac_height = 72*2.54*cm;
  
  G4ThreeVector dewar_inner_vac_pos = G4ThreeVector(0*cm, 0*cm, -dewar_inner_vac_height/2.0 + dewar_top_to_MC + he_det_center_to_MC);
  G4Tubs* dewar_inner_vac_shape =    
     new G4Tubs("dewar_inner_vac", 
         dewar_inner_vac_OD/2.0 - dewar_inner_vac_thick, 
         dewar_inner_vac_OD/2.0, 
         dewar_inner_vac_height/2.0, 
         0.*deg, 
         360.*deg);
  
  G4LogicalVolume* dewar_inner_vac_logic = 
     new G4LogicalVolume(dewar_inner_vac_shape,
                         ss_t316,
                         "dewar_inner_vac");
  
  dewar_inner_vac_logic->SetVisAttributes(redTVA);
  
  new G4PVPlacement(0,                       //no rotation
                    dewar_inner_vac_pos,     //at position
                    dewar_inner_vac_logic,   //its logical volume
                    "dewar_inner_vac_",      //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  //dewar helium bath
  G4double dewar_bath_helium_OD = dewar_inner_vac_OD - 2.0*dewar_inner_vac_thick;
  G4double dewar_bath_helium_ID = ivc_can_OD;
  G4double dewar_bath_helium_height = ivc_can_height;
  
  G4ThreeVector dewar_bath_helium_pos = ivc_can_pos;
  G4Tubs* dewar_bath_helium_shape =    
     new G4Tubs("dewar_bath_helium", 
         dewar_bath_helium_ID/2.0, 
         dewar_bath_helium_OD/2.0, 
         dewar_bath_helium_height/2.0, 
         0.*deg, 
         360.*deg);
  
  G4LogicalVolume* dewar_bath_helium_logic = 
     new G4LogicalVolume(dewar_bath_helium_shape,
                         liquid_helium,
                         "dewar_bath_helium");
  dewar_bath_helium_logic->SetVisAttributes(helium_vis);
  
  new G4PVPlacement(0,                       //no rotation
                    dewar_bath_helium_pos,//at position
                    dewar_bath_helium_logic,   //its logical volume
                    "dewar_bath_helium",      //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking  
  
  //dewar inner nitrogen can
  G4double dewar_inner_nitrogen_OD = 18.062*2.54*cm;
  G4double dewar_inner_nitrogen_thick = 0.105*2.54*cm;
  G4double dewar_inner_nitrogen_height = 60.5*2.54*cm;
  
  G4ThreeVector dewar_inner_nitrogen_pos = G4ThreeVector(0*cm, 0*cm, -9.5*2.54*cm - dewar_inner_nitrogen_height/2.0 + dewar_top_to_MC + he_det_center_to_MC);
  G4Tubs* dewar_inner_nitrogen_shape =    
     new G4Tubs("dewar_inner_nitrogen", 
         dewar_inner_nitrogen_OD/2.0 - dewar_inner_nitrogen_thick, 
         dewar_inner_nitrogen_OD/2.0, 
         dewar_inner_nitrogen_height/2.0, 
         0.*deg, 
         360.*deg);
  
  G4LogicalVolume* dewar_inner_nitrogen_logic = 
     new G4LogicalVolume(dewar_inner_nitrogen_shape,
                         ss_t316,
                         "dewar_inner_nitrogen");
  
  dewar_inner_nitrogen_logic->SetVisAttributes(redTVA);
  
  new G4PVPlacement(0,                       //no rotation
                    dewar_inner_nitrogen_pos,//at position
                    dewar_inner_nitrogen_logic,   //its logical volume
                    "dewar_inner_nitrogen",      //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking  

  //dewar outer nitrogen can
  G4double dewar_outer_nitrogen_OD = 20.105*2.54*cm;
  G4double dewar_outer_nitrogen_thick = 0.105*2.54*cm;
  G4double dewar_outer_nitrogen_height = 60.5*2.54*cm;
  
  G4ThreeVector dewar_outer_nitrogen_pos = G4ThreeVector(0*cm, 0*cm, -9.5*2.54*cm - dewar_outer_nitrogen_height/2.0 + dewar_top_to_MC + he_det_center_to_MC);
  G4Tubs* dewar_outer_nitrogen_shape =    
     new G4Tubs("dewar_outer_nitrogen", 
         dewar_outer_nitrogen_OD/2.0 - dewar_inner_nitrogen_thick, 
         dewar_outer_nitrogen_OD/2.0, 
         dewar_outer_nitrogen_height/2.0, 
         0.*deg, 
         360.*deg);
  
  G4LogicalVolume* dewar_outer_nitrogen_logic = 
     new G4LogicalVolume(dewar_outer_nitrogen_shape,
                         ss_t316,
                         "dewar_outer_nitrogen");
  
  dewar_outer_nitrogen_logic->SetVisAttributes(redTVA);
  
  new G4PVPlacement(0,                       //no rotation
                    dewar_outer_nitrogen_pos,//at position
                    dewar_outer_nitrogen_logic,   //its logical volume
                    "dewar_outer_nitrogen",      //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking  

  //dewar nitrogen bath
  G4double dewar_bath_nitrogen_OD = dewar_outer_nitrogen_OD - 2.0*dewar_outer_nitrogen_thick;
  G4double dewar_bath_nitrogen_ID = dewar_inner_nitrogen_OD;
  G4double dewar_bath_nitrogen_height = dewar_inner_nitrogen_height;
  
  G4ThreeVector dewar_bath_nitrogen_pos = dewar_inner_nitrogen_pos;
  G4Tubs* dewar_bath_nitrogen_shape =    
     new G4Tubs("dewar_bath_nitrogen", 
         dewar_bath_nitrogen_ID/2.0, 
         dewar_bath_nitrogen_OD/2.0, 
         dewar_bath_nitrogen_height/2.0, 
         0.*deg, 
         360.*deg);
  
  G4LogicalVolume* dewar_bath_nitrogen_logic = 
     new G4LogicalVolume(dewar_bath_nitrogen_shape,
                         liquid_nitrogen,
                         "dewar_bath_nitrogen");
  
  dewar_bath_nitrogen_logic->SetVisAttributes(blueTVA);
  
  new G4PVPlacement(0,                       //no rotation
                    dewar_bath_nitrogen_pos, //at position
                    dewar_bath_nitrogen_logic,//its logical volume
                    "dewar_bath_nitrogen",   //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking  

  //dewar outer vacuum can
  G4double dewar_outer_vacuum_OD = 22.0*2.54*cm;
  G4double dewar_outer_vacuum_thick = 0.105*2.54*cm;
  G4double dewar_outer_vacuum_height = 87.375*2.54*cm;
  
  G4ThreeVector dewar_outer_vacuum_pos = G4ThreeVector(0*cm, 0*cm, -5.0*2.54*cm - dewar_outer_vacuum_height/2.0 + dewar_top_to_MC + he_det_center_to_MC);
  G4Tubs* dewar_outer_vacuum_shape =    
     new G4Tubs("dewar_outer_vacuum", 
         dewar_outer_vacuum_OD/2.0 - dewar_outer_vacuum_thick, 
         dewar_outer_vacuum_OD/2.0, 
         dewar_outer_vacuum_height/2.0, 
         0.*deg, 
         360.*deg);
  
  G4LogicalVolume* dewar_outer_vacuum_logic = 
     new G4LogicalVolume(dewar_outer_vacuum_shape,
                         ss_t316,
                         "dewar_outer_vacuum");
  
  dewar_outer_vacuum_logic->SetVisAttributes(redTVA);
  
  new G4PVPlacement(0,                       //no rotation
                    dewar_outer_vacuum_pos, //at position
                    dewar_outer_vacuum_logic,//its logical volume
                    "dewar_outer_vacuum",    //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking  


  std::cout << "Detector Created" << std::endl;
  
  //
  // Always return the physical World
  //
  return physWorld;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void B2aDetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors

  G4String trackerChamberSDname = "B2/TrackerChamberSD";
  B2TrackerSD* aTrackerSD = new B2TrackerSD(trackerChamberSDname,
                                            "TrackerHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);
  // Setting aTrackerSD to all logical volumes with the same name 
  // of "Chamber_LV".
  //SetSensitiveDetector("Chamber_LV", aTrackerSD, true);

  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue = G4ThreeVector();
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void B2aDetectorConstruction::SetTargetMaterial(G4String materialName)
{
  G4NistManager* nistManager = G4NistManager::Instance();

  G4Material* pttoMaterial = 
              nistManager->FindOrBuildMaterial(materialName);

  if (fTargetMaterial != pttoMaterial) {
     if ( pttoMaterial ) {
        fTargetMaterial = pttoMaterial;
        if (fLogicTarget) fLogicTarget->SetMaterial(fTargetMaterial);
        G4cout 
          << G4endl 
          << "----> The target is made of " << materialName << G4endl;
     } else {
        G4cout 
          << G4endl 
          << "-->  WARNING from SetTargetMaterial : "
          << materialName << " not found" << G4endl;
     }
  }
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2aDetectorConstruction::SetChamberMaterial(G4String materialName)
{
  G4NistManager* nistManager = G4NistManager::Instance();

  G4Material* pttoMaterial =
              nistManager->FindOrBuildMaterial(materialName);

  if (fChamberMaterial != pttoMaterial) {
     if ( pttoMaterial ) {
        fChamberMaterial = pttoMaterial;
        for (G4int copyNo=0; copyNo<fNbOfChambers; copyNo++) {
            if (fLogicChamber[copyNo]) fLogicChamber[copyNo]->
                                               SetMaterial(fChamberMaterial);
        }
        G4cout 
          << G4endl 
          << "----> The chambers are made of " << materialName << G4endl;
     } else {
        G4cout 
          << G4endl 
          << "-->  WARNING from SetChamberMaterial : "
          << materialName << " not found" << G4endl;
     }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2aDetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((fStepLimit)&&(maxStep>0.)) fStepLimit->SetMaxAllowedStep(maxStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2aDetectorConstruction::SetCheckOverlaps(G4bool checkOverlaps)
{
  fCheckOverlaps = checkOverlaps;
}  
