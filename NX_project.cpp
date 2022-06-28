#include <stdio.h>
#include <uf.h> 
#include <uf_csys.h> 
#include <uf_curve.h> 
#include <uf_defs.h>
#include <uf_mtx.h> 
#include <uf_ui.h>
#include <cmath>
#include <NXOpen/Annotations_Dimension.hxx>
#include <NXOpen/Arc.hxx>
#include <NXOpen/Axis.hxx>
#include <NXOpen/AxisCollection.hxx>
#include <NXOpen/BasePart.hxx>
#include <NXOpen/Body.hxx>
#include <NXOpen/BodyCollection.hxx>
#include <NXOpen/BodyDumbRule.hxx>
#include <NXOpen/Builder.hxx>
#include <NXOpen/CartesianCoordinateSystem.hxx>
#include <NXOpen/CoordinateSystem.hxx>
#include <NXOpen/CoordinateSystemCollection.hxx>
#include <NXOpen/CurveCollection.hxx>
#include <NXOpen/CurveDumbRule.hxx>
#include <NXOpen/CurveFeatureRule.hxx>
#include <NXOpen/DatumAxis.hxx>
#include <NXOpen/DatumCollection.hxx>
#include <NXOpen/DatumPlane.hxx>
#include <NXOpen/Direction.hxx>
#include <NXOpen/DirectionCollection.hxx>
#include <NXOpen/DisplayableObject.hxx>
#include <NXOpen/DisplayManager.hxx>
#include <NXOpen/DisplayModification.hxx>
#include <NXOpen/Edge.hxx>
#include <NXOpen/Expression.hxx>
#include <NXOpen/ExpressionCollection.hxx>
#include <NXOpen/Face.hxx>
#include <NXOpen/FaceTangentRule.hxx>
#include <NXOpen/Features_BooleanBuilder.hxx>
#include <NXOpen/Features_Cone.hxx>
#include <NXOpen/Features_ConeBuilder.hxx>
#include <NXOpen/Features_CylinderBuilder.hxx>
#include <NXOpen/Features_DatumCsys.hxx>
#include <NXOpen/Features_DatumCsysBuilder.hxx>
#include <NXOpen/Features_Feature.hxx>
#include <NXOpen/Features_FeatureCollection.hxx>
#include <NXOpen/Features_Helix.hxx>
#include <NXOpen/Features_HelixBuilder.hxx>
#include <NXOpen/Features_LawCurve.hxx>
#include <NXOpen/Features_LawCurveBuilder.hxx>
#include <NXOpen/Features_SketchFeature.hxx>
#include <NXOpen/Features_SphereBuilder.hxx>
#include <NXOpen/Features_Swept.hxx>
#include <NXOpen/Features_SweptBuilder.hxx>
#include <NXOpen/GeometricUtilities_AlignmentMethodBuilder.hxx>
#include <NXOpen/GeometricUtilities_AlongSpineBuilder.hxx>
#include <NXOpen/GeometricUtilities_BooleanRegionSelect.hxx>
#include <NXOpen/GeometricUtilities_FeatureOptions.hxx>
#include <NXOpen/GeometricUtilities_LawBuilder.hxx>
#include <NXOpen/GeometricUtilities_MultiTransitionLawBuilder.hxx>
#include <NXOpen/GeometricUtilities_NonInflectingLawBuilder.hxx>
#include <NXOpen/GeometricUtilities_OnPathDimensionBuilder.hxx>
#include <NXOpen/GeometricUtilities_OrientationMethodBuilder.hxx>
#include <NXOpen/GeometricUtilities_Rebuild.hxx>
#include <NXOpen/GeometricUtilities_ScalingMethodBuilder.hxx>
#include <NXOpen/GeometricUtilities_SShapedLawBuilder.hxx>
#include <NXOpen/IBaseCurve.hxx>
#include <NXOpen/IPlane.hxx>
#include <NXOpen/IReferenceAxis.hxx>
#include <NXOpen/Line.hxx>
#include <NXOpen/ModelingView.hxx>
#include <NXOpen/ModelingViewCollection.hxx>
#include <NXOpen/NXException.hxx>
#include <NXOpen/NXMatrix.hxx>
#include <NXOpen/NXObject.hxx>
#include <NXOpen/ObjectList.hxx>
#include <NXOpen/Offset.hxx>
#include <NXOpen/OffsetCollection.hxx>
#include <NXOpen/Part.hxx>
#include <NXOpen/PartCollection.hxx>
#include <NXOpen/Plane.hxx>
#include <NXOpen/PlaneCollection.hxx>
#include <NXOpen/PlaneTypes.hxx>
#include <NXOpen/Point.hxx>
#include <NXOpen/PointCollection.hxx>
#include <NXOpen/Preferences_RulePreferences.hxx>
#include <NXOpen/Preferences_SessionPreferences.hxx>
#include <NXOpen/Preferences_SessionSketch.hxx>
#include <NXOpen/Preferences_SketchPreferences.hxx>
#include <NXOpen/PreviewBuilder.hxx>
#include <NXOpen/RuleManager.hxx>
#include <NXOpen/Scalar.hxx>
#include <NXOpen/ScalarCollection.hxx>
#include <NXOpen/ScCollector.hxx>
#include <NXOpen/ScCollectorCollection.hxx>
#include <NXOpen/ScRuleFactory.hxx>
#include <NXOpen/Section.hxx>
#include <NXOpen/SectionCollection.hxx>
#include <NXOpen/SectionList.hxx>
#include <NXOpen/SelectBodyList.hxx>
#include <NXOpen/SelectDisplayableObjectList.hxx>
#include <NXOpen/SelectFaceList.hxx>
#include <NXOpen/SelectionIntentRule.hxx>
#include <NXOpen/SelectIReferenceAxis.hxx>
#include <NXOpen/SelectISurface.hxx>
#include <NXOpen/SelectLine.hxx>
#include <NXOpen/SelectObject.hxx>
#include <NXOpen/SelectObjectList.hxx>
#include <NXOpen/Session.hxx>
#include <NXOpen/Sketch.hxx>
#include <NXOpen/SketchAlongPathBuilder.hxx>
#include <NXOpen/SketchCollection.hxx>
#include <NXOpen/SketchDimensionalConstraint.hxx>
#include <NXOpen/SketchGeometricConstraint.hxx>
#include <NXOpen/SketchInPlaceBuilder.hxx>
#include <NXOpen/SmartObject.hxx>
#include <NXOpen/Spline.hxx>
#include <NXOpen/TaggedObject.hxx>
#include <NXOpen/Unit.hxx>
#include <NXOpen/UnitCollection.hxx>
#include <NXOpen/Update.hxx>
#include <NXOpen/View.hxx>
#include <NXOpen/WCS.hxx>
#include <NXOpen/Xform.hxx>
#include <NXOpen/XformCollection.hxx>


NXOpen::Point3d RKGetLineIntersectionPoint(NXOpen::Point3d& a, NXOpen::Point3d& b, NXOpen::Point3d& c, NXOpen::Point3d& d);
void ufusr(char* param, int* retcode, int paramLen)
{
    
   // double d1 = {47.68};
  //  double d2;
    double d3 = {47.68};
    double ln = {76};
    double phi = { 7.125*PI/180 };
    double phi_grad = { 7.125};
    double P = { 5.08 };
    double h1 = {2.993};
    double f = {0.508};
    double r = { 0.508 };
    double r1 = { 0.38 };
	
	if (UF_initialize()) return;
	
	
    NXOpen::Session* theSession = NXOpen::Session::GetSession();
    NXOpen::Part* workPart(theSession->Parts()->Work());
    NXOpen::Part* displayPart(theSession->Parts()->Display());
    
    //конус---------------------------------------------------------------------------------------------------------------

    
    NXOpen::Features::Cone* nullNXOpen_Features_Cone(NULL);
    NXOpen::Features::ConeBuilder* coneBuilder1;
    coneBuilder1 = workPart->Features()->CreateConeBuilder(nullNXOpen_Features_Cone);

    coneBuilder1->SetType(NXOpen::Features::ConeBuilder::TypesTopDiameterHeightAndHalfAngle);

    coneBuilder1->TopDiameter()->SetValue(d3-0.0052);

    coneBuilder1->Height()->SetValue(ln);

    coneBuilder1->HalfAngle()->SetValue(phi_grad);

    coneBuilder1->BooleanOption()->SetType(NXOpen::GeometricUtilities::BooleanOperation::BooleanTypeCreate);

      
    NXOpen::NXObject* cone1;
    cone1 = coneBuilder1->Commit();	 
   
  
    coneBuilder1->Destroy();
    	    
    NXOpen::Features::Cone* con1(dynamic_cast<NXOpen::Features::Cone*>(cone1));
    std::vector <NXOpen::Face*> FaceAr = con1->GetFaces();

	//цилиндр************************************************************
       

    NXOpen::Features::Feature* nullNXOpen_Features_Feature(NULL);
    NXOpen::Features::CylinderBuilder* cylinderBuilder1;
    cylinderBuilder1 = workPart->Features()->CreateCylinderBuilder(nullNXOpen_Features_Feature);       
    cylinderBuilder1->Diameter()->SetFormula("80");
    cylinderBuilder1->Height()->SetFormula("80");
    
    //объединение цилиндра с конусом---------------------------------------------------------------------------------------------------------------
    
    cylinderBuilder1->BooleanOption()->SetType(NXOpen::GeometricUtilities::BooleanOperation::BooleanTypeUnite);
    
    std::vector<NXOpen::Body*> targetBodies2(1);
    NXOpen::Body* body1(FaceAr[0]->GetBody());// = dynamic_cast<NXOpen::Body*>(c);
    
    targetBodies2[0] = body1;
    cylinderBuilder1->BooleanOption()->SetTargetBodies(targetBodies2);
   
    
    NXOpen::DatumAxis* datumAxis1(dynamic_cast<NXOpen::DatumAxis*>(workPart->Datums()->FindObject("DATUM_CSYS(0) Z axis")));
    NXOpen::Direction* direction1;
    direction1 = workPart->Directions()->CreateDirection(datumAxis1, NXOpen::SenseForward, NXOpen::SmartObject::UpdateOptionWithinModeling);

    NXOpen::Point* nullNXOpen_Point(NULL);
    NXOpen::Axis* axis1;
    axis1 = workPart->Axes()->CreateAxis(nullNXOpen_Point, direction1, NXOpen::SmartObject::UpdateOptionWithinModeling);

    cylinderBuilder1->SetAxis(axis1);
    axis1->SetPoint(nullNXOpen_Point);
   // axis1->Evaluate();
    direction1->ReverseDirection();
    axis1->SetDirection(direction1);                   
    NXOpen::Point3d p1(0.0, 0.0, 0.0);
    NXOpen::Point* point4;
    point4 = workPart->Points()->CreatePoint(p1);    
    axis1->SetPoint(point4);         
    NXOpen::NXObject* cyl1;
    cyl1 = cylinderBuilder1->Commit();      
    cylinderBuilder1->Destroy();
     
    
    //эскиз закона изменения диаметра спирали---------------------------------------------------------------------------------------------------------------
    
    NXOpen::Point3d origin1(0.0, 0.0, 0.0);
    NXOpen::Point3d axis2(0.0, 1.0, 0.0);
    NXOpen::Matrix3x3 wcs_matrix = workPart->WCS()->CoordinateSystem()->Orientation()->Element();
    NXOpen::DatumPlane* datumPlane2 = workPart->Datums()->CreateFixedDatumPlane(origin1, wcs_matrix);
    NXOpen::DatumAxis* datumAxis2 = workPart->Datums()->CreateFixedDatumAxis(origin1, axis2);

    NXOpen::Sketch* nullNXOpen_Sketch(NULL);
    NXOpen::SketchInPlaceBuilder* sketchInPlaceBuilder1;
    sketchInPlaceBuilder1 = workPart->Sketches()->CreateSketchInPlaceBuilder2(nullNXOpen_Sketch);


    sketchInPlaceBuilder1->PlaneOrFace()->SetValue(datumPlane2);
    sketchInPlaceBuilder1->Axis()->SetValue(datumAxis2);
    sketchInPlaceBuilder1->SetSketchOrigin(workPart->Points()->CreatePoint(origin1));
    sketchInPlaceBuilder1->SetPlaneOption(NXOpen::Sketch::PlaneOption::PlaneOptionInferred);

    NXOpen::NXObject* ob1;
    ob1 = sketchInPlaceBuilder1->Commit();

    sketchInPlaceBuilder1->Destroy();

    NXOpen::Sketch* sk1(dynamic_cast<NXOpen::Sketch*>(ob1));
    sk1->Activate(NXOpen::Sketch::ViewReorientTrue);
   
    //базовая линия
    NXOpen::Point3d pl1(0, 0, 0);
    NXOpen::Point3d p2(50, 0, 0);
    NXOpen::Line *l1= workPart->Curves()->CreateLine(pl1, p2);
    theSession->ActiveSketch()->AddGeometry(l1, NXOpen::Sketch::InferConstraintsOptionInferNoConstraints);
    //основные линии
    NXOpen::Point3d p3(0, d3/2, 0);
    NXOpen::Point3d p4(ln-1.5*P, d3/2+ (ln - 1.5 * P)*sin(phi), 0);
    NXOpen::Point3d p5(ln, d3/2 + (ln - 1.5 * P) * sin(phi)+1.5*P*sin(PI/6), 0);
    NXOpen::Line* l2 = workPart->Curves()->CreateLine(p3, p4);
    NXOpen::Line* l3 = workPart->Curves()->CreateLine(p4, p5);
    theSession->ActiveSketch()->AddGeometry(l2, NXOpen::Sketch::InferConstraintsOptionInferNoConstraints);
    theSession->ActiveSketch()->AddGeometry(l3, NXOpen::Sketch::InferConstraintsOptionInferNoConstraints);
    theSession->ActiveSketch()->Deactivate(NXOpen::Sketch::ViewReorientTrue, NXOpen::Sketch::UpdateLevelModel);
    


    //система координат для спирали---------------------------------------------------------------------------------------------------------------------
       
   
    NXOpen::Features::Feature* nullNXOpen_Features_Feature1(NULL);
    NXOpen::Features::DatumCsysBuilder* datumCsysBuilder1;
    datumCsysBuilder1 = workPart->Features()->CreateDatumCsysBuilder(nullNXOpen_Features_Feature1);

    

    NXOpen::Unit* unit1(dynamic_cast<NXOpen::Unit*>(workPart->UnitCollection()->FindObject("MilliMeter")));
    NXOpen::Unit* unit2(dynamic_cast<NXOpen::Unit*>(workPart->UnitCollection()->FindObject("Degrees")));
    NXOpen::Expression* expression1;
    expression1 = workPart->Expressions()->CreateSystemExpressionWithUnits("0", unit1);

  

    NXOpen::Point3d origin4(0.0, 0.0, 0.0);
    NXOpen::Vector3d normal4(0.0, 0.0, 1.0);
    NXOpen::Plane* plane4;
    plane4 = workPart->Planes()->CreatePlane(origin4, normal4, NXOpen::SmartObject::UpdateOptionWithinModeling);

  

    NXOpen::Expression* expression12;
    expression12 = workPart->Expressions()->CreateSystemExpressionWithUnits("0", unit1);

    NXOpen::Expression* expression13;
    expression13 = workPart->Expressions()->CreateSystemExpressionWithUnits("0", unit2);

    NXOpen::Expression* expression14;
    expression14 = workPart->Expressions()->CreateSystemExpressionWithUnits("0", unit1);

    NXOpen::Expression* expression15;
    expression15 = workPart->Expressions()->CreateSystemExpressionWithUnits("0", unit2);

    NXOpen::Expression* expression16;
    expression16 = workPart->Expressions()->CreateSystemExpressionWithUnits("0", unit1);

    NXOpen::Expression* expression17;
    expression17 = workPart->Expressions()->CreateSystemExpressionWithUnits("0", unit2);

   
    expression12->SetFormula("0");

    expression14->SetFormula("0");

    expression16->SetFormula("76");

    expression13->SetFormula("180");

    expression15->SetFormula("0");

    expression17->SetFormula("0");

    expression12->SetRightHandSide("0");

    expression14->SetRightHandSide("0");

    expression16->SetRightHandSide("76");

    expression13->SetRightHandSide("180");

    expression15->SetRightHandSide("0");

    expression17->SetRightHandSide("0");   

    expression12->SetFormula("0");

    expression14->SetFormula("0");

    expression16->SetFormula("76");

    expression13->SetFormula("180");

    expression15->SetFormula("0");

    expression17->SetFormula("0");

  
    NXOpen::Point3d origin5(0.0, 0.0, 0.0);
    NXOpen::Matrix3x3 orientation1;
    orientation1.Xx = 1.0;
    orientation1.Xy = 0.0;
    orientation1.Xz = 0.0;
    orientation1.Yx = 0.0;
    orientation1.Yy = 1.0;
    orientation1.Yz = 0.0;
    orientation1.Zx = 0.0;
    orientation1.Zy = 0.0;
    orientation1.Zz = 1.0;
    NXOpen::CartesianCoordinateSystem* cartesianCoordinateSystem1;
    cartesianCoordinateSystem1 = workPart->CoordinateSystems()->CreateCoordinateSystem(origin5, orientation1, true);

    NXOpen::Point3d origin6(0.0, 0.0, 0.0);
    NXOpen::Matrix3x3 orientation2;
    orientation2.Xx = 1.0;
    orientation2.Xy = 0.0;
    orientation2.Xz = 0.0;
    orientation2.Yx = 0.0;
    orientation2.Yy = 1.0;
    orientation2.Yz = 0.0;
    orientation2.Zx = 0.0;
    orientation2.Zy = 0.0;
    orientation2.Zz = 1.0;
    NXOpen::CartesianCoordinateSystem* cartesianCoordinateSystem2;
    cartesianCoordinateSystem2 = workPart->CoordinateSystems()->CreateCoordinateSystem(origin6, orientation2, true);

    NXOpen::Scalar* scalar1;
    scalar1 = workPart->Scalars()->CreateScalarExpression(expression12, NXOpen::Scalar::DimensionalityTypeLength, NXOpen::SmartObject::UpdateOptionWithinModeling);

    NXOpen::Scalar* scalar2;
    scalar2 = workPart->Scalars()->CreateScalarExpression(expression14, NXOpen::Scalar::DimensionalityTypeLength, NXOpen::SmartObject::UpdateOptionWithinModeling);

    NXOpen::Scalar* scalar3;
    scalar3 = workPart->Scalars()->CreateScalarExpression(expression16, NXOpen::Scalar::DimensionalityTypeLength, NXOpen::SmartObject::UpdateOptionWithinModeling);

    NXOpen::Scalar* scalar4;
    scalar4 = workPart->Scalars()->CreateScalarExpression(expression13, NXOpen::Scalar::DimensionalityTypeAngle, NXOpen::SmartObject::UpdateOptionWithinModeling);

    NXOpen::Scalar* scalar5;
    scalar5 = workPart->Scalars()->CreateScalarExpression(expression15, NXOpen::Scalar::DimensionalityTypeAngle, NXOpen::SmartObject::UpdateOptionWithinModeling);

    NXOpen::Scalar* scalar6;
    scalar6 = workPart->Scalars()->CreateScalarExpression(expression17, NXOpen::Scalar::DimensionalityTypeAngle, NXOpen::SmartObject::UpdateOptionWithinModeling);

    NXOpen::Offset* offset1;
    offset1 = workPart->Offsets()->CreateOffsetRectangular(scalar1, scalar2, scalar3, NXOpen::SmartObject::UpdateOptionWithinModeling);

    NXOpen::Vector3d originOffset1(0.0, 0.0, 0.0);
    NXOpen::Matrix3x3 trasformMatrix1;
    trasformMatrix1.Xx = 1.0;
    trasformMatrix1.Xy = 0.0;
    trasformMatrix1.Xz = 0.0;
    trasformMatrix1.Yx = 0.0;
    trasformMatrix1.Yy = 1.0;
    trasformMatrix1.Yz = 0.0;
    trasformMatrix1.Zx = 0.0;
    trasformMatrix1.Zy = 0.0;
    trasformMatrix1.Zz = 1.0;
    NXOpen::Xform* xform1;
    xform1 = workPart->Xforms()->CreateXformByDynamicOffset(cartesianCoordinateSystem2, originOffset1, trasformMatrix1, NXOpen::SmartObject::UpdateOptionWithinModeling, 1.0);

    NXOpen::CartesianCoordinateSystem* cartesianCoordinateSystem3;
    cartesianCoordinateSystem3 = workPart->CoordinateSystems()->CreateCoordinateSystem(xform1, NXOpen::SmartObject::UpdateOptionWithinModeling);

    NXOpen::Offset* nullNXOpen_Offset(NULL);
    NXOpen::Xform* xform2;
    xform2 = workPart->Xforms()->CreateXform(cartesianCoordinateSystem3, nullNXOpen_Offset, offset1, scalar4, scalar5, scalar6, 0, NXOpen::SmartObject::UpdateOptionWithinModeling, 1.0);

   
   
    NXOpen::CartesianCoordinateSystem* cartesianCoordinateSystem4;
    cartesianCoordinateSystem4 = workPart->CoordinateSystems()->CreateCoordinateSystem(xform2, NXOpen::SmartObject::UpdateOptionWithinModeling);  



    plane4->DestroyPlane();

    datumCsysBuilder1->SetCsys(cartesianCoordinateSystem4);

    datumCsysBuilder1->SetDisplayScaleFactor(1.25);

    NXOpen::NXObject* nXObject1;
    nXObject1 = datumCsysBuilder1->Commit();

    datumCsysBuilder1->Destroy();

    
    //спираль по закону, заданному прямыми линиями---------------------------------------------------------------------------------------------------------------
    
   
    NXOpen::Features::Helix* nullNXOpen_Features_Helix(NULL);
    NXOpen::Features::HelixBuilder* helixBuilder1;
    helixBuilder1 = workPart->Features()->CreateHelixBuilder(nullNXOpen_Features_Helix);
   

    helixBuilder1->SetOrientationOption(NXOpen::Features::HelixBuilder::OrientationOptionsSpecified);

    helixBuilder1->StartAngle()->SetFormula("180");

    helixBuilder1->SetSizeOption(NXOpen::Features::HelixBuilder::SizeOptionsRadius);

    helixBuilder1->SizeLaw()->SetLawType(NXOpen::GeometricUtilities::LawBuilder::TypeByLawCurve);

    
    helixBuilder1->PitchLaw()->Value()->SetFormula("5.08");

    helixBuilder1->PitchLaw()->StartValue()->SetFormula("5");

    helixBuilder1->PitchLaw()->EndValue()->SetFormula("5");

   

    helixBuilder1->StartLimit()->SetPercentUsed(false);

    helixBuilder1->StartLimit()->Expression()->SetFormula("0");

    helixBuilder1->EndLimit()->SetPercentUsed(false);

    helixBuilder1->EndLimit()->Expression()->SetFormula("76");
    

    

   /* NXOpen::DatumPlane* datumPlane1(dynamic_cast<NXOpen::DatumPlane*>(nXObject1));
    NXOpen::Xform* xform10;
    xform10 = workPart->Xforms()->CreateXform(nXObject1, NXOpen::SmartObject::UpdateOptionWithinModeling);

    NXOpen::CartesianCoordinateSystem* cartesianCoordinateSystem10;
    cartesianCoordinateSystem10 = workPart->CoordinateSystems()->CreateCoordinateSystem(xform10, NXOpen::SmartObject::UpdateOptionWithinModeling);
    */
    
    helixBuilder1->SetCoordinateSystem(cartesianCoordinateSystem4);
    
    helixBuilder1->SizeLaw()->LawCurve()->SetAllowedEntityTypes(NXOpen::Section::AllowTypesOnlyCurves);

    
    
    std::vector<NXOpen::IBaseCurve*> curves1(1);
    NXOpen::Line* line1 = l2;
   // curves1[0] = l5;  
    NXOpen::CurveDumbRule* curveDumbRule1;
    curveDumbRule1 = workPart->ScRuleFactory()->CreateRuleBaseCurveDumb({ l2 });
    

    helixBuilder1->SizeLaw()->LawCurve()->AllowSelfIntersection(true);

    std::vector<NXOpen::SelectionIntentRule*> rules1(1);
    rules1[0] = curveDumbRule1;
    NXOpen::NXObject* nullNXOpen_NXObject(NULL);
    NXOpen::Point3d helpPoint1(0,0,0);
    helixBuilder1->SizeLaw()->LawCurve()->AddToSection(rules1, l2, nullNXOpen_NXObject, nullNXOpen_NXObject, helpPoint1, NXOpen::Section::ModeCreate, false);

    
   
  
    helixBuilder1->Evaluate();

  
   
    std::vector<NXOpen::IBaseCurve*> curves2(1);
    NXOpen::Line* line2(l3);
    curves2[0] = line2;
    NXOpen::CurveDumbRule* curveDumbRule2;
    curveDumbRule2 = workPart->ScRuleFactory()->CreateRuleBaseCurveDumb({ l3 });

    helixBuilder1->SizeLaw()->LawCurve()->AllowSelfIntersection(true);
    
    std::vector<NXOpen::SelectionIntentRule *> rules2(1);
    rules2[0] = curveDumbRule2;
    NXOpen::Point3d helpPoint2(0,0,0);
    helixBuilder1->SizeLaw()->LawCurve()->AddToSection(rules2, l3, nullNXOpen_NXObject, nullNXOpen_NXObject, helpPoint2, NXOpen::Section::ModeCreate, false);
  
  
   
    NXOpen::Line* line3(l1);
    helixBuilder1->SizeLaw()->BaseLine()->SetValue(line3);

    helixBuilder1->SizeLaw()->SetReverseDirection(true);

    helixBuilder1->Evaluate();

    helixBuilder1->Evaluate();

    helixBuilder1->SizeLaw()->SetReverseDirection(false);

    helixBuilder1->Evaluate();

    helixBuilder1->EndLimit()->Expression()->SetFormula("75");

    helixBuilder1->EndLimit()->Expression()->SetFormula("75");

    
    NXOpen::NXObject* hel1;
    hel1 = helixBuilder1->Commit();
     
     
    helixBuilder1->Destroy();
    
    
 
    //эскиз профиля резьбы--------------------------------------------------------------------------------------------------------
    NXOpen::Point3d origin2(0.0, 0.0, 0.0);
    NXOpen::Point3d axis3(0.0, 0.0, 1.0);
    NXOpen::Matrix3x3 wcs_matrix1 = workPart->WCS()->CoordinateSystem()->Orientation()->Element();
    double mt[9] = { wcs_matrix1.Xx,wcs_matrix1.Xy,wcs_matrix1.Xz, wcs_matrix1.Yx,wcs_matrix1.Yy,wcs_matrix1.Yz, wcs_matrix1.Zx,wcs_matrix1.Zy,wcs_matrix1.Zz };
    double mtP[9];
    double vec[3] = { 1,0,0 };
    UF_MTX3_rotate_about_axis(vec, 90 * PI / 180,mtP);
    UF_MTX3_multiply(mtP, mt, mt);
    wcs_matrix1.Xx = mt[0];
    wcs_matrix1.Xy = mt[1];
    wcs_matrix1.Xz = mt[2];
    wcs_matrix1.Yx = mt[3];
    wcs_matrix1.Yy = mt[4]; 
    wcs_matrix1.Yz = mt[5]; 
    wcs_matrix1.Zx = mt[6]; 
    wcs_matrix1.Zy = mt[7];
    wcs_matrix1.Zz = mt[8];
    NXOpen::DatumPlane* datumPlane3 = workPart->Datums()->CreateFixedDatumPlane(origin2, wcs_matrix1);
    NXOpen::DatumAxis* datumAxis3 = workPart->Datums()->CreateFixedDatumAxis(origin2, axis3);


    NXOpen::SketchInPlaceBuilder* sketchInPlaceBuilder2;
    sketchInPlaceBuilder2 = workPart->Sketches()->CreateSketchInPlaceBuilder2(NULL);


    sketchInPlaceBuilder2->PlaneOrFace()->SetValue(datumPlane3);
    sketchInPlaceBuilder2->Axis()->SetValue(datumAxis3);
    sketchInPlaceBuilder2->SetSketchOrigin(workPart->Points()->CreatePoint(origin2));
    sketchInPlaceBuilder2->SetPlaneOption(NXOpen::Sketch::PlaneOption::PlaneOptionInferred);

    NXOpen::NXObject* ob2;
    ob2 = sketchInPlaceBuilder2->Commit();

    sketchInPlaceBuilder2->Destroy();

    NXOpen::Sketch* sk2(dynamic_cast<NXOpen::Sketch*>(ob2));
    sk2->Activate(NXOpen::Sketch::ViewReorientTrue);

    
    NXOpen::Point3d tp1(-1 * d3 / 2, 0 , ln);
    NXOpen::Point3d tp9(-1 * d3 / 2 +P/2*sin(phi)+h1+f, 0 , ln +P/2);
    NXOpen::Point3d tp6(-1 * d3 / 2 + P * sin(phi), 0, ln + P);
    NXOpen::Point3d tp7(-1 * d3 / 2 + P * sin(phi)-1, 0, ln + P);
    NXOpen::Point3d tp8(-1 * d3 / 2-1, 0, ln);


    NXOpen::Point3d temp1(-1 * d3 / 2 + P / 2 * sin(phi) + h1 + f-r * cos(30 * DEGRA) /tan(30*PI/180), 0, ln + P / 2-r*cos(30*DEGRA));
    NXOpen::Point3d temp2(-1 * d3 / 2 + P / 2 * sin(phi) + h1 + f - r * cos(30 * DEGRA) / tan(30 * PI / 180), 0, ln + P / 2 + r * cos(30 * DEGRA));
    NXOpen::Point3d tp2 = RKGetLineIntersectionPoint(tp1, tp6, tp9, temp1);
    NXOpen::Point3d tp5 = RKGetLineIntersectionPoint(tp1, tp6, tp9, temp2);
    
    double h3 = 0.38 / (sin((120 + phi_grad)/2*PI/180));
    NXOpen::Point3d tc1(tp2.X+h3*sin((60- phi_grad /2)*PI/180), tp2.Y, tp2.Z- h3 * cos((60 - phi_grad / 2) * PI / 180));
    NXOpen::Point3d tf1(tc1.X-0.38*cos(phi), tc1.Y, tc1.Z+0.38*sin(phi));
    NXOpen::Point3d tf2(tc1.X - 0.38 * sin(PI/6), tc1.Y, tc1.Z + 0.38 * cos(PI/6));
    NXOpen::Point3d tf3(tc1.X - 0.38 * sin(PI / 4), tc1.Y, tc1.Z + 0.38 * cos(PI / 4));
    
    double h4 = 0.38 / (sin((120 - phi_grad) / 2 * PI / 180));
    NXOpen::Point3d tc2(tp5.X + h4 * sin((60 + phi_grad / 2) * PI / 180), tp5.Y, tp5.Z + h4 * cos((60 + phi_grad / 2) * PI / 180));
    NXOpen::Point3d tf4(tc2.X - 0.38 * cos(phi), tc2.Y, tc2.Z + 0.38 * sin(phi));
    NXOpen::Point3d tf5(tc2.X - 0.38 * sin(PI / 6), tc2.Y, tc2.Z - 0.38 * cos(PI / 6));
    NXOpen::Point3d tf6(tc2.X - 0.38 * sin(PI / 4), tc2.Y, tc2.Z - 0.38 * cos(PI / 4));

    NXOpen::Point3d tc3(temp1.X-r*sin(30*DEGRA)+r, 0, ln + P / 2);
  
    bool b1 = { FALSE };
    bool* bool1;
    bool1 = &b1;
    NXOpen::Arc* arc1;
    arc1 = workPart->Curves()->CreateArc(tf1, tf3, tf2, false,bool1);
    
    bool b2 = { FALSE };
    bool* bool2;
    bool2 = &b2;
    NXOpen::Arc* arc2;
    arc2 = workPart->Curves()->CreateArc(tf4, tf6, tf5, false, bool2);

    bool b3 = { FALSE };
    bool* bool3;
    bool3 = &b3;
    NXOpen::Arc* arc3;
    arc3 = workPart->Curves()->CreateArc(temp1, tc3, temp2, false, bool3);


    theSession->ActiveSketch()->AddGeometry(arc1, NXOpen::Sketch::InferConstraintsOptionInferNoConstraints);
    theSession->ActiveSketch()->AddGeometry(arc2, NXOpen::Sketch::InferConstraintsOptionInferNoConstraints);
    theSession->ActiveSketch()->AddGeometry(arc3, NXOpen::Sketch::InferConstraintsOptionInferNoConstraints);
    
    NXOpen::Line* tl1 = workPart->Curves()->CreateLine(tp1, tf1);
    NXOpen::Line* tl2 = workPart->Curves()->CreateLine(tf2, temp1);
   NXOpen::Line* tl3 = workPart->Curves()->CreateLine(temp2, tf5);
   NXOpen::Line* tl4 = workPart->Curves()->CreateLine(tf4, tp6);
   NXOpen::Line* tl5 = workPart->Curves()->CreateLine(tp6, tp7);
   NXOpen::Line* tl6 = workPart->Curves()->CreateLine(tp7, tp8);
   NXOpen::Line* tl7 = workPart->Curves()->CreateLine(tp8, tp1);

    theSession->ActiveSketch()->AddGeometry(tl1, NXOpen::Sketch::InferConstraintsOptionInferNoConstraints);
    theSession->ActiveSketch()->AddGeometry(tl2, NXOpen::Sketch::InferConstraintsOptionInferNoConstraints);
    theSession->ActiveSketch()->AddGeometry(tl3, NXOpen::Sketch::InferConstraintsOptionInferNoConstraints);
    theSession->ActiveSketch()->AddGeometry(tl4, NXOpen::Sketch::InferConstraintsOptionInferNoConstraints);
    theSession->ActiveSketch()->AddGeometry(tl5, NXOpen::Sketch::InferConstraintsOptionInferNoConstraints);
    theSession->ActiveSketch()->AddGeometry(tl6, NXOpen::Sketch::InferConstraintsOptionInferNoConstraints);
    theSession->ActiveSketch()->AddGeometry(tl7, NXOpen::Sketch::InferConstraintsOptionInferNoConstraints);
  
       theSession->ActiveSketch()->Deactivate(NXOpen::Sketch::ViewReorientTrue, NXOpen::Sketch::UpdateLevelModel);

   
  //заметание
       
       NXOpen::Features::Swept* nullNXOpen_Features_Swept(NULL);
       NXOpen::Features::SweptBuilder* sweptBuilder1;
       sweptBuilder1 = workPart->Features()->CreateSweptBuilder(nullNXOpen_Features_Swept);

       sweptBuilder1->SetG0Tolerance(0.01);

       sweptBuilder1->SetG1Tolerance(0.5);

       sweptBuilder1->OrientationMethod()->SetOrientationOption(NXOpen::GeometricUtilities::OrientationMethodBuilder::OrientationOptionsByFaceNormals);

       sweptBuilder1->OrientationMethod()->AngularLaw()->Value()->SetFormula("7.125");

       sweptBuilder1->OrientationMethod()->AngularLaw()->StartValue()->SetFormula("0");

       sweptBuilder1->OrientationMethod()->AngularLaw()->EndValue()->SetFormula("0");

       sweptBuilder1->ScalingMethod()->AreaLaw()->Value()->SetFormula("1");

       sweptBuilder1->ScalingMethod()->AreaLaw()->StartValue()->SetFormula("1");

       sweptBuilder1->ScalingMethod()->AreaLaw()->EndValue()->SetFormula("1");

       sweptBuilder1->ScalingMethod()->PerimeterLaw()->Value()->SetFormula("1");

       sweptBuilder1->ScalingMethod()->PerimeterLaw()->StartValue()->SetFormula("1");

       sweptBuilder1->ScalingMethod()->PerimeterLaw()->EndValue()->SetFormula("1");
       
       sweptBuilder1->Spine()->SetDistanceTolerance(0.01);

       sweptBuilder1->Spine()->SetChainingTolerance(0.0094999999999999998);

       sweptBuilder1->AlignmentMethod()->AlignCurve()->SetDistanceTolerance(0.01);

       sweptBuilder1->AlignmentMethod()->AlignCurve()->SetChainingTolerance(0.0094999999999999998);

       sweptBuilder1->OrientationMethod()->OrientationCurve()->SetDistanceTolerance(0.01);

       sweptBuilder1->OrientationMethod()->OrientationCurve()->SetChainingTolerance(0.0094999999999999998);

       sweptBuilder1->OrientationMethod()->AngularLaw()->AlongSpineData()->Spine()->SetDistanceTolerance(0.01);

       sweptBuilder1->OrientationMethod()->AngularLaw()->AlongSpineData()->Spine()->SetChainingTolerance(0.0094999999999999998);

       sweptBuilder1->OrientationMethod()->AngularLaw()->LawCurve()->SetDistanceTolerance(0.01);

       sweptBuilder1->OrientationMethod()->AngularLaw()->LawCurve()->SetChainingTolerance(0.0094999999999999998);

       sweptBuilder1->ScalingMethod()->ScalingCurve()->SetDistanceTolerance(0.01);

       sweptBuilder1->ScalingMethod()->ScalingCurve()->SetChainingTolerance(0.0094999999999999998);

       sweptBuilder1->ScalingMethod()->AreaLaw()->AlongSpineData()->Spine()->SetDistanceTolerance(0.01);

       sweptBuilder1->ScalingMethod()->AreaLaw()->AlongSpineData()->Spine()->SetChainingTolerance(0.0094999999999999998);

       sweptBuilder1->ScalingMethod()->AreaLaw()->LawCurve()->SetDistanceTolerance(0.01);

       sweptBuilder1->ScalingMethod()->AreaLaw()->LawCurve()->SetChainingTolerance(0.0094999999999999998);

       sweptBuilder1->ScalingMethod()->PerimeterLaw()->AlongSpineData()->Spine()->SetDistanceTolerance(0.01);

       sweptBuilder1->ScalingMethod()->PerimeterLaw()->AlongSpineData()->Spine()->SetChainingTolerance(0.0094999999999999998);

       sweptBuilder1->ScalingMethod()->PerimeterLaw()->LawCurve()->SetDistanceTolerance(0.01);

       sweptBuilder1->ScalingMethod()->PerimeterLaw()->LawCurve()->SetChainingTolerance(0.0094999999999999998);

       sweptBuilder1->Spine()->SetAngleTolerance(0.5);

       sweptBuilder1->AlignmentMethod()->AlignCurve()->SetAngleTolerance(0.5);

       sweptBuilder1->OrientationMethod()->OrientationCurve()->SetAngleTolerance(0.5);

       sweptBuilder1->OrientationMethod()->AngularLaw()->AlongSpineData()->Spine()->SetAngleTolerance(0.5);

       sweptBuilder1->OrientationMethod()->AngularLaw()->LawCurve()->SetAngleTolerance(0.5);

       sweptBuilder1->ScalingMethod()->ScalingCurve()->SetAngleTolerance(0.5);

       sweptBuilder1->ScalingMethod()->AreaLaw()->AlongSpineData()->Spine()->SetAngleTolerance(0.5);

       sweptBuilder1->ScalingMethod()->AreaLaw()->LawCurve()->SetAngleTolerance(0.5);

       sweptBuilder1->ScalingMethod()->PerimeterLaw()->AlongSpineData()->Spine()->SetAngleTolerance(0.5);

       sweptBuilder1->ScalingMethod()->PerimeterLaw()->LawCurve()->SetAngleTolerance(0.5);

       NXOpen::Section* section1;
       section1 = workPart->Sections()->CreateSection(0.0094999999999999998, 0.01, 0.5);

       sweptBuilder1->SectionList()->Append(section1);

       section1->SetAllowedEntityTypes(NXOpen::Section::AllowTypesOnlyCurves);

       std::vector<NXOpen::Features::Feature*> features1(1);
       //NXOpen::Features::SketchFeature* sketchFeature1(dynamic_cast<NXOpen::Features::SketchFeature*>(ob2));
     //  NXOpen::Features::SketchFeature* sketchFeature1(ob2);
      // features1[0] = sketchFeature1;
       NXOpen::CurveDumbRule* curveFeatureRule1;
       curveFeatureRule1 = workPart->ScRuleFactory()->CreateRuleCurveDumb({ tl1,arc1, tl2,arc3, tl3, arc2, tl4, tl5, tl6, tl7});
       
       section1->AllowSelfIntersection(false);

       std::vector<NXOpen::SelectionIntentRule*> rules3(1);
       rules3[0] = curveFeatureRule1;
      // NXOpen::NXObject* nullNXOpen_NXObject(NULL);
       //NXOpen::Point3d helpPoint1(0.0, 0.0, 0.0);
       section1->AddToSection(rules3, nullNXOpen_NXObject, nullNXOpen_NXObject, nullNXOpen_NXObject, helpPoint1, NXOpen::Section::ModeCreate, false);

       std::vector<NXOpen::Section*> sections1(1);
       sections1[0] = section1;
       sweptBuilder1->AlignmentMethod()->SetSections(sections1);

       NXOpen::Point3d pointonstartcurve1(-24.34, 0.0, 76.0);
       section1->SetStartCurveOfClosedLoop(0, pointonstartcurve1);

       section1->ReverseDirectionOfClosedLoop(0);

       try
       {
           sweptBuilder1->AlignmentMethod()->UpdateSectionAtIndex(0);
       }
       catch (const NXOpen::NXException& ex)
       {
           ex.AssertErrorCode(1);
       }
       
       NXOpen::Section* section2;
       section2 = workPart->Sections()->CreateSection(0.0094999999999999998, 0.01, 0.5);

       sweptBuilder1->GuideList()->Append(section2);

       section2->SetAllowedEntityTypes(NXOpen::Section::AllowTypesOnlyCurves);

       std::vector<NXOpen::Features::Feature*> features2(1);
       NXOpen::Features::Helix* helix1(dynamic_cast<NXOpen::Features::Helix*>(hel1));
       features2[0] = helix1;
       NXOpen::CurveFeatureRule* curveFeatureRule2;
       curveFeatureRule2 = workPart->ScRuleFactory()->CreateRuleCurveFeature(features2);
       
       section2->AllowSelfIntersection(false);

       std::vector<NXOpen::SelectionIntentRule*> rules4(1);
       rules4[0] = curveFeatureRule2;
      // NXOpen::Point3d helpPoint2(0.0, 0.0, 0.0);
       section2->AddToSection(rules4, nullNXOpen_NXObject, nullNXOpen_NXObject, nullNXOpen_NXObject, helpPoint2, NXOpen::Section::ModeCreate, false);

       sweptBuilder1->ScalingMethod()->AreaLaw()->AlongSpineData()->SetFeatureSpine(section2);

       sweptBuilder1->ScalingMethod()->PerimeterLaw()->AlongSpineData()->SetFeatureSpine(section2);

       sweptBuilder1->OrientationMethod()->AngularLaw()->AlongSpineData()->SetFeatureSpine(section2);
      // NXOpen::Body* bd(dynamic_cast<NXOpen::Body* > (cone1));
      
      
       
      // NXOpen::Face* face1(dynamic_cast<NXOpen::Face*>(con1->FindObject("FACE 2 {(-28.5899889875762,-0,38) CONE(1)}")));
       NXOpen::Face* face1(FaceAr[0]);
       std::vector<NXOpen::Face*> boundaryFaces1(0);
       NXOpen::FaceTangentRule* faceTangentRule1;
       faceTangentRule1 = workPart->ScRuleFactory()->CreateRuleFaceTangent(face1, boundaryFaces1, 0.5);

       std::vector<NXOpen::SelectionIntentRule*> rules5(1);
       rules5[0] = faceTangentRule1;
       sweptBuilder1->OrientationMethod()->Faces()->ReplaceRules(rules5, false);

       NXOpen::NXObject* sw1;
       sw1 = sweptBuilder1->Commit();
       
       sweptBuilder1->Destroy();
       

       //вычитание конуса и заметания


       NXOpen::Features::BooleanFeature* nullNXOpen_Features_BooleanFeature(NULL);
       NXOpen::Features::BooleanBuilder* booleanBuilder1;
       booleanBuilder1 = workPart->Features()->CreateBooleanBuilderUsingCollector(nullNXOpen_Features_BooleanFeature);

       NXOpen::ScCollector* scCollector1;
       scCollector1 = booleanBuilder1->ToolBodyCollector();

       NXOpen::GeometricUtilities::BooleanRegionSelect* booleanRegionSelect1;
       booleanRegionSelect1 = booleanBuilder1->BooleanRegionSelect();

       booleanBuilder1->SetTolerance(0.01);

       booleanBuilder1->SetOperation(NXOpen::Features::Feature::BooleanTypeSubtract);

       
       //NXOpen::Body* sub_body1(dynamic_cast<NXOpen::Body*>(cone1));
       NXOpen::Body* sub_body1 = FaceAr[0]->GetBody();
       bool added1;
       added1 = booleanBuilder1->Targets()->Add(sub_body1);
       
       
       std::vector<NXOpen::TaggedObject*> targets1(1);
       targets1[0] = sub_body1;
       booleanRegionSelect1->AssignTargets(targets1);
       
       NXOpen::ScCollector* scCollector2;
       scCollector2 = workPart->ScCollectors()->CreateCollector();

       NXOpen::Features::Swept* swep1(dynamic_cast<NXOpen::Features::Swept*>(sw1));
       std::vector <NXOpen::Face*> FaceAr2 = swep1->GetFaces();

       std::vector<NXOpen::Body*> bodies1(1);
       NXOpen::Body* sub_body2(FaceAr2[0]->GetBody());
       bodies1[0] = sub_body2;
       NXOpen::BodyDumbRule* bodyDumbRule1;
       bodyDumbRule1 = workPart->ScRuleFactory()->CreateRuleBodyDumb(bodies1, true);
       
       std::vector<NXOpen::SelectionIntentRule*> sub_rules1(1);
       sub_rules1[0] = bodyDumbRule1;
       scCollector2->ReplaceRules(sub_rules1, false);

       booleanBuilder1->SetToolBodyCollector(scCollector2);

       std::vector<NXOpen::TaggedObject*> targets2(1);
       targets2[0] = sub_body1;
       booleanRegionSelect1->AssignTargets(targets2);

            
      
       NXOpen::NXObject* sub1;
      sub1 = booleanBuilder1->Commit();
     
      
       booleanBuilder1->Destroy();

       //скрытие элементов дерева построения

       std::vector<NXOpen::DisplayableObject*> dis_objects2(6);     
       dis_objects2[0] = sk1;
       dis_objects2[1] = sk2;
       dis_objects2[2] = datumPlane3;
       dis_objects2[3] = datumAxis3;
       dis_objects2[4] = datumPlane2;
       dis_objects2[5] = datumAxis2;
       theSession->DisplayManager()->BlankObjects(dis_objects2);


	UF_terminate();
}
int ufusr_ask_unload(void) 
{
	return (UF_UNLOAD_IMMEDIATELY);
}

NXOpen::Point3d RKGetLineIntersectionPoint(NXOpen::Point3d& a, NXOpen::Point3d& b, NXOpen::Point3d& c, NXOpen::Point3d& d)
{
    double k1 = (b.X - a.X) / (b.Z - a.Z);
    double m1 = -1 * a.Z * k1 + a.X;
    double k2 = (d.X - c.X) / (d.Z - c.Z);
    double m2 = -1 * c.Z * k2 + c.X;
    NXOpen::Point3d res;
    res.Z = (m1 - m2) / (k2 - k1);
    res.X = k1 * (m1 - m2) / (k2 - k1) + m1;
    res.Y = 0;
    return res;
}
