#include "Joint.h"
Joint::Joint(NXOpen::Session* theSession, NXOpen::Part* workPart, NXOpen::Part* displayPart) {
    d1 = { 47.68 };
    d3 = { 41.0 };
    ln = { 60 };
    phi = { 5.7105 * DEGRA };
    phi_grad = { 5.7105 };
    P = { 4.233 };
    h1 = { 2.5 };
    f = { 0.423 };
    f1 = { 0.731 };
    r = { 0.508 };
    r1 = { 0.38 };
    h = { 2.192 };
    d4 = { 54 };
    d5 = { 48.616 };
    pin_cyl_diam = { 64 };
    pin_cyl_len = { 20 };
    box_cyl_len = { 90 };
    thr_esc = { 12.7 };
    box_l1 = { 67 };
    lm = { 78 };
    cone_recess_len = { 16 };//конусна€ выточка в муфте, длина
    inner_hole_diam = { 22 };
    anglerot = 10;
    rotZ_Offset = -anglerot / 360 * P / 2;
    fric_coef = 0.13;
    diam_avg_thread = d3 + 2 * tan(phi) * (ln - 12.7);//диметр резьбы в главном сечении
    diam_avg_endFace = (d3 + 2 * tan(phi) * (ln)+pin_cyl_diam) / 2; //средний диаметр торца
    cfactor << std::scientific << ln / 86 * 1.752E-08 * 0;
    mesh_size = std::to_string(1.3);
    Joint::theSession= theSession;
    Joint::workPart = workPart;
    Joint::displayPart = displayPart;
   
}
void Joint::build(){
    NXOpen::NXObject* cone1 = buildCone(d3 - 0.005, ln, phi_grad, { 0,0,0 }, { 0,0,1 }, 0, nullptr, 0);

    NXOpen::Features::Cone* con1(dynamic_cast<NXOpen::Features::Cone*>(cone1));
    FaceAr = con1->GetFaces();
    
    //цилиндр************************************************************


    NXOpen::NXObject* cyl1 = buildCylinder(pin_cyl_diam, pin_cyl_len, { 0,0,0 }, { 0,0,-1 }, cone1, 1);


    //эскиз закона изменени€ диаметра спирали---------------------------------------------------------------------------------------------------------------
    NXOpen::Sketch* sk1(dynamic_cast<NXOpen::Sketch*>(buildSketch({ 0,0,0 }, { 0,0,1 }, { 1,0,0 })));
    sk1->Activate(NXOpen::Sketch::ViewReorientTrue);


    //базова€ лини€
    NXOpen::Point3d pl1(0, 0, 0);
    NXOpen::Point3d p2(50, 0, 0);
    NXOpen::Line* l1 = workPart->Curves()->CreateLine(pl1, p2);
    theSession->ActiveSketch()->AddGeometry(l1, NXOpen::Sketch::InferConstraintsOptionInferNoConstraints);
    //основные линии ниппел€
    NXOpen::Point3d p3(0, d3 / 2, 0);
    NXOpen::Point3d p4(ln - thr_esc, d3 / 2 + (ln - thr_esc) * tan(phi), 0);
    NXOpen::Point3d p5(ln - thr_esc + 1.5 * P, d3 / 2 + (ln - thr_esc) * tan(phi) + 1.5 * P * tan(PI / 6), 0);
    NXOpen::Line* l2 = workPart->Curves()->CreateLine(p3, p4);
    NXOpen::Line* l3 = workPart->Curves()->CreateLine(p4, p5);

    //основные линии муфты
    NXOpen::Point3d hel2p1(0, d3 / 2, 0);
    NXOpen::Point3d hel2p2(-box_l1, d3 / 2 + box_l1 * tan(phi), 0);

    NXOpen::Point3d hel2p3(1.5 * P, d3 / 2 - 1.5 * P * tan(30 * DEGRA), 0);

    NXOpen::Line* hel2l1 = workPart->Curves()->CreateLine(hel2p1, hel2p2);
    NXOpen::Line* hel2l2 = workPart->Curves()->CreateLine(hel2p1, hel2p3);

    theSession->ActiveSketch()->AddGeometry(l2, NXOpen::Sketch::InferConstraintsOptionInferNoConstraints);
    theSession->ActiveSketch()->AddGeometry(l3, NXOpen::Sketch::InferConstraintsOptionInferNoConstraints);
    theSession->ActiveSketch()->AddGeometry(hel2l1, NXOpen::Sketch::InferConstraintsOptionInferNoConstraints);
    theSession->ActiveSketch()->AddGeometry(hel2l2, NXOpen::Sketch::InferConstraintsOptionInferNoConstraints);



    theSession->ActiveSketch()->Deactivate(NXOpen::Sketch::ViewReorientTrue, NXOpen::Sketch::UpdateLevelModel);



    //система координат дл€ спирали---------------------------------------------------------------------------------------------------------------------


    NXOpen::CartesianCoordinateSystem* cartesianCoordinateSystem4 = buildCsys({ 0,0,ln }, { 1,0,0 }, { 0,-1,0 });
    //спираль по закону, заданному пр€мыми лини€ми---------------------------------------------------------------------------------------------------------------

    NXOpen::Line* hel_lines1[] = { l2,l3 };
    NXOpen::NXObject* hel1 = helixByLawCurve(90.0, P, 0.0, ln - thr_esc + 1.5 * P, cartesianCoordinateSystem4, l1, hel_lines1, 2, false);

    NXOpen::CartesianCoordinateSystem* cartesianCoordinateSystem5 = buildCsys({ 0,0,box_l1 + rotZ_Offset }, { 1,0,0 }, { 0,-1,0 });
    NXOpen::Line* hel_lines2[] = { hel2l1,hel2l2 };
    NXOpen::NXObject* hel2 = helixByLawCurve(270, P, -1.5 * P, box_l1, cartesianCoordinateSystem5, l1, hel_lines2, 2, true);


    //эскиз профил€ резьбы--------------------------------------------------------------------------------------------------------




    NXOpen::Point3d p1[8];

    p1[0] = { 0,0,0 };
    p1[1] = { 0,f1,0 };
    p1[6] = { P,P * sin(phi) + f1,0 };
    p1[7] = { P,P * sin(phi),0 };
    p1[2] = getLinesIntersectionPoint(p1[0], 60, p1[1], phi_grad);
    p1[3] = getLinesIntersectionPoint({ 0,f1 + h1,0 }, phi_grad, p1[0], 60);
    p1[4] = getLinesIntersectionPoint({ 0,f1 + h1,0 }, phi_grad, p1[7], 120);
    p1[5] = getLinesIntersectionPoint(p1[7], 120, p1[1], phi_grad);


    NXOpen::Point3d points2[8];
    for (int i = 0; i < 8; i++) {
        points2[i] = p1[i];

    }

    NXOpen::Point3d points3[8];
    NXOpen::Point3d ttp2 = getLinesIntersectionPoint(p1[0], 60, p1[7], 120);
    points3[0] = multiplyByMatrix4x4(p1[0], { -1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1 });
    double offsetX = ttp2.X - points3[0].X;
    double offsetY = ttp2.Y - points3[0].Y;
    double offsetZ = ttp2.Z - points3[0].Z;
    for (int i = 0; i < 8; i++) {
        points3[i] = multiplyByMatrix4x4(p1[i], { -1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1 });
        points3[i] = multiplyByMatrix4x4(points3[i], { 1,0,0,offsetX,0,1,0,offsetY,0,0,1,offsetZ,0,0,0,1 });
    }


    for (int i = 0; i < 8; i++) {
        points2[i] = multiplyByMatrix4x4(points2[i], { 1,0,0,ln,0,1,0,-d3 / 2 - f1,0,0,1,0,0,0,0,1 });
        points2[i] = multiplyByMatrix4x4(points2[i], { 0,0,-1,0,0,1,0,0,1,0,0,0,0,0,0,1 });
        points3[i] = multiplyByMatrix4x4(points3[i], { 1,0,0,ln,0,1,0,-d3 / 2 - f1,0,0,1,0,0,0,0,1 });
        points3[i] = multiplyByMatrix4x4(points3[i], { 0,0,-1,0,0,1,0,0,1,0,0,0,0,0,0,1 });
        points3[i] = multiplyByMatrix4x4(points3[i], { cos(anglerot * DEGRA),sin(anglerot * DEGRA),0,0,-sin(anglerot * DEGRA),cos(anglerot * DEGRA),0,0,0,0,1,0,0,0,0,1 });
        points3[i] = multiplyByMatrix4x4(points3[i], { 1,0,0,0,0,1,0,0,0,0,1,rotZ_Offset,0,0,0,1 });
    }



    NXOpen::Sketch* sk2(dynamic_cast<NXOpen::Sketch*>(buildSketch({ 0,0,0 }, { -1,0,0 }, { 0,0,1 })));
    sk2->Activate(NXOpen::Sketch::ViewReorientTrue);

    NXOpen::Line* lines1[8];
    lines1[0] = workPart->Curves()->CreateLine(points2[0], points2[1]);
    lines1[1] = workPart->Curves()->CreateLine(points2[1], points2[2]);
    lines1[2] = workPart->Curves()->CreateLine(points2[2], points2[3]);
    lines1[3] = workPart->Curves()->CreateLine(points2[3], points2[4]);
    lines1[4] = workPart->Curves()->CreateLine(points2[4], points2[5]);
    lines1[5] = workPart->Curves()->CreateLine(points2[5], points2[6]);
    lines1[6] = workPart->Curves()->CreateLine(points2[6], points2[7]);
    lines1[7] = workPart->Curves()->CreateLine(points2[7], points2[0]);
    for (int i = 0; i < 8; i++) {
        theSession->ActiveSketch()->AddGeometry(lines1[i], NXOpen::Sketch::InferConstraintsOptionInferNoConstraints);
    }


    NXOpen::Line* lines2[8];
    lines2[0] = workPart->Curves()->CreateLine(points3[0], points3[1]);
    lines2[1] = workPart->Curves()->CreateLine(points3[1], points3[2]);
    lines2[2] = workPart->Curves()->CreateLine(points3[2], points3[3]);
    lines2[3] = workPart->Curves()->CreateLine(points3[3], points3[4]);
    lines2[4] = workPart->Curves()->CreateLine(points3[4], points3[5]);
    lines2[5] = workPart->Curves()->CreateLine(points3[5], points3[6]);
    lines2[6] = workPart->Curves()->CreateLine(points3[6], points3[7]);
    lines2[7] = workPart->Curves()->CreateLine(points3[7], points3[0]);

    theSession->ActiveSketch()->Deactivate(NXOpen::Sketch::ViewReorientTrue, NXOpen::Sketch::UpdateLevelModel);

    NXOpen::Sketch* sk3(dynamic_cast<NXOpen::Sketch*>(buildSketch({ 0,0,0 }, { -1,0,0 }, { 0,0,1 })));

    sk3->Activate(NXOpen::Sketch::ViewReorientTrue);
    for (int i = 0; i < 8; i++) {
        theSession->ActiveSketch()->AddGeometry(lines2[i], NXOpen::Sketch::InferConstraintsOptionInferNoConstraints);
    }
    theSession->ActiveSketch()->Deactivate(NXOpen::Sketch::ViewReorientTrue, NXOpen::Sketch::UpdateLevelModel);


    //заметание
    NXOpen::NXObject* sw1 = nullptr;
    bool swept_flag = true;
    while (swept_flag) {

        NXOpen::Features::Swept* nullNXOpen_Features_Swept(NULL);
        NXOpen::Features::SweptBuilder* sweptBuilder1;
        sweptBuilder1 = workPart->Features()->CreateSweptBuilder(nullNXOpen_Features_Swept);

        sweptBuilder1->SetG0Tolerance(0.01);
        sweptBuilder1->SetG1Tolerance(0.5);
        sweptBuilder1->OrientationMethod()->SetOrientationOption(NXOpen::GeometricUtilities::OrientationMethodBuilder::OrientationOptionsByFaceNormals);
        /*sweptBuilder1->OrientationMethod()->AngularLaw()->Value()->SetFormula("7.125");*/
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
        NXOpen::CurveDumbRule* curveFeatureRule1;
        curveFeatureRule1 = workPart->ScRuleFactory()->CreateRuleCurveDumb({ lines1[0],lines1[1],lines1[2],lines1[3],lines1[4],lines1[5],lines1[6],lines1[7] });

        section1->AllowSelfIntersection(false);

        std::vector<NXOpen::SelectionIntentRule*> rules3(1);
        rules3[0] = curveFeatureRule1;
        NXOpen::NXObject* nullNXOpen_NXObject(NULL);
        NXOpen::Point3d helpPoint1(0.0, 0.0, 0.0);
        section1->AddToSection(rules3, nullNXOpen_NXObject, nullNXOpen_NXObject, nullNXOpen_NXObject, helpPoint1, NXOpen::Section::ModeCreate, false);

        std::vector<NXOpen::Section*> sections1(1);
        sections1[0] = section1;
        sweptBuilder1->AlignmentMethod()->SetSections(sections1);

        NXOpen::Point3d pointonstartcurve1(points2[0]);
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
        section2->AddToSection(rules4, nullNXOpen_NXObject, nullNXOpen_NXObject, nullNXOpen_NXObject, { 0,0,0 }, NXOpen::Section::ModeCreate, false);

        sweptBuilder1->ScalingMethod()->AreaLaw()->AlongSpineData()->SetFeatureSpine(section2);
        sweptBuilder1->ScalingMethod()->PerimeterLaw()->AlongSpineData()->SetFeatureSpine(section2);
        sweptBuilder1->OrientationMethod()->AngularLaw()->AlongSpineData()->SetFeatureSpine(section2);

        NXOpen::Face* face1(FaceAr[0]);
        std::vector<NXOpen::Face*> boundaryFaces1(0);
        NXOpen::FaceTangentRule* faceTangentRule1;
        faceTangentRule1 = workPart->ScRuleFactory()->CreateRuleFaceTangent(face1, boundaryFaces1, 0.5);

        std::vector<NXOpen::SelectionIntentRule*> rules5(1);
        rules5[0] = faceTangentRule1;
        sweptBuilder1->OrientationMethod()->Faces()->ReplaceRules(rules5, false);
        _sleep(1000);
        try {
            sw1 = sweptBuilder1->Commit();
        }
        catch (...) { sweptBuilder1->Destroy(); UF_terminate(); };
        sweptBuilder1->Destroy();
        swept_flag = false;

    }

    //вычитание конуса и заметани€

    NXOpen::NXObject* ob5 = subtractBodies(cone1, sw1);

    //cyl2
    NXOpen::NXObject* cyl2 = buildCylinder(pin_cyl_diam, box_cyl_len, { 0,0, rotZ_Offset }, { 0,0,1 }, nullptr, 0);

    //cone2
    NXOpen::NXObject* cone2 = buildCone(d5 + 0.03, lm, phi_grad, { 0,0, rotZ_Offset }, { 0,0,1 }, 1, cyl2, 2);

    NXOpen::Features::Cone* con2(dynamic_cast<NXOpen::Features::Cone*>(cone2));
    FaceAr1 = con2->GetFaces();

    NXOpen::NXObject* sw2;
    {

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
        NXOpen::CurveDumbRule* curveFeatureRule1;
        curveFeatureRule1 = workPart->ScRuleFactory()->CreateRuleCurveDumb({ lines2[0],lines2[1],lines2[2],lines2[3],lines2[4],lines2[5],lines2[6],lines2[7] });

        section1->AllowSelfIntersection(false);

        std::vector<NXOpen::SelectionIntentRule*> rules3(1);
        rules3[0] = curveFeatureRule1;
        NXOpen::NXObject* nullNXOpen_NXObject(NULL);
        NXOpen::Point3d helpPoint1(0.0, 0.0, 0.0);
        section1->AddToSection(rules3, nullNXOpen_NXObject, nullNXOpen_NXObject, nullNXOpen_NXObject, helpPoint1, NXOpen::Section::ModeCreate, false);

        std::vector<NXOpen::Section*> sections1(1);
        sections1[0] = section1;
        sweptBuilder1->AlignmentMethod()->SetSections(sections1);

        NXOpen::Point3d pointonstartcurve1(points3[0]);
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
        NXOpen::Features::Helix* helix1(dynamic_cast<NXOpen::Features::Helix*>(hel2));
        features2[0] = helix1;
        NXOpen::CurveFeatureRule* curveFeatureRule2;
        curveFeatureRule2 = workPart->ScRuleFactory()->CreateRuleCurveFeature(features2);

        section2->AllowSelfIntersection(false);

        std::vector<NXOpen::SelectionIntentRule*> rules4(1);
        rules4[0] = curveFeatureRule2;
        // NXOpen::Point3d helpPoint2(0.0, 0.0, 0.0);
        section2->AddToSection(rules4, nullNXOpen_NXObject, nullNXOpen_NXObject, nullNXOpen_NXObject, { 0,0,0 }, NXOpen::Section::ModeCreate, false);

        sweptBuilder1->ScalingMethod()->AreaLaw()->AlongSpineData()->SetFeatureSpine(section2);
        sweptBuilder1->ScalingMethod()->PerimeterLaw()->AlongSpineData()->SetFeatureSpine(section2);
        sweptBuilder1->OrientationMethod()->AngularLaw()->AlongSpineData()->SetFeatureSpine(section2);

        NXOpen::Face* face1(FaceAr1[0]);
        std::vector<NXOpen::Face*> boundaryFaces1(0);
        NXOpen::FaceTangentRule* faceTangentRule1;
        faceTangentRule1 = workPart->ScRuleFactory()->CreateRuleFaceTangent(face1, boundaryFaces1, 0.5);

        std::vector<NXOpen::SelectionIntentRule*> rules5(1);
        rules5[0] = faceTangentRule1;
        sweptBuilder1->OrientationMethod()->Faces()->ReplaceRules(rules5, false);


        sw2 = sweptBuilder1->Commit();

        sweptBuilder1->Destroy();
    }


    NXOpen::NXObject* ob6 = subtractBodies(cone2, sw2);

    NXOpen::NXObject* cone3 = buildCone(d4, cone_recess_len, phi_grad, { 0,0,rotZ_Offset }, { 0,0,1 }, 1, ob6, 2);
    NXOpen::NXObject* cone4 = buildCone((d4 / 2 - cone_recess_len * tan(phi)) * 2, cone_recess_len, 30, { 0,0,cone_recess_len + rotZ_Offset }, { 0,0,1 }, 1, cone3, 2);
    NXOpen::NXObject* cyl3 = buildCylinder(inner_hole_diam, box_cyl_len * 2, { 0,0,box_cyl_len + rotZ_Offset }, { 0,0,-1 }, cone4, 2);
    NXOpen::NXObject* cyl4 = buildCylinder(inner_hole_diam, box_cyl_len * 2, { 0,0,box_cyl_len + rotZ_Offset }, { 0,0,-1 }, cone1, 2);


    int type_face;
    int count = 0;
    double parms[2], face_point[3];
    double point_on_face[3] = { 0,0,0 };

    NXOpen::Features::BodyFeature* b1(dynamic_cast<NXOpen::Features::BodyFeature*>(cyl1));
    std::vector <NXOpen::Face*> cylFacesAr1 = b1->GetFaces();
    for (auto& face : cylFacesAr1)
    {
        UF_MODL_ask_face_type(face->Tag(), &type_face);
        if (type_face == UF_MODL_CYLINDRICAL_FACE)
        {
            face->SetColor(108);//цилинлрическа€ поверхность ниппел€
            count++;
        }
        if (type_face == UF_MODL_PLANAR_FACE)
        {
            UF_MODL_ask_face_parm_2(face->Tag(), point_on_face, parms, face_point);
            if (face_point[2] >= -1) face->SetColor(6);//плоска€ поверхность ниппел€
        }


    }


    point_on_face[2] = { 200 };
    NXOpen::Features::BodyFeature* b2(dynamic_cast<NXOpen::Features::BodyFeature*>(cyl2));
    std::vector <NXOpen::Face*> cylFacesAr2 = b2->GetFaces();
    for (auto& face : cylFacesAr2)
    {
        UF_MODL_ask_face_type(face->Tag(), &type_face);
        UF_MODL_ask_face_parm_2(face->Tag(), point_on_face, parms, face_point);
        if (type_face == UF_MODL_PLANAR_FACE)
        {
            UF_MODL_ask_face_parm_2(face->Tag(), point_on_face, parms, face_point);
            if (face_point[2] >= 5) face->SetColor(211);//дальн€€ плоска€ поверхность муты
            else face->SetColor(83);//ближн€€ плоска€ поверхность муфты

            count++;
        }

    }


    NXOpen::Features::BodyFeature* b3(dynamic_cast<NXOpen::Features::BodyFeature*>(sw1));
    std::vector <NXOpen::Face*> cylFacesAr3 = b3->GetFaces();
    for (auto& face : cylFacesAr3)
    {
       
        face->SetColor(36);
        

    }


    NXOpen::Features::BodyFeature* b4(dynamic_cast<NXOpen::Features::BodyFeature*>(sw2));
    std::vector <NXOpen::Face*> cylFacesAr4 = b4->GetFaces();
    for (auto& face : cylFacesAr4)
    {
      
        face->SetColor(83);
        

    }


    //скрытие элементов дерева построени€

    std::vector<NXOpen::DisplayableObject*> dis_objects2(6);
    dis_objects2[0] = sk1;
    dis_objects2[1] = sk2;
    dis_objects2[2] = sk3;
    dis_objects2[3] = dynamic_cast<NXOpen::Spline*> (hel1);
    dis_objects2[4] = dynamic_cast<NXOpen::Spline*> (hel2);


    theSession->DisplayManager()->BlankObjects(dis_objects2);
}

/**
Creates a cylinder
params:
origin - start point
direction - direction vector
body2 - target body to boolean operation
boolType (0 - Create without boolean operation
    1 - Unite objects
    2 - Subtract objects
    3 - Intersect objects)
*/
NXOpen::NXObject* Joint::buildCylinder(double diameter, double height, NXOpen::Point3d origin, NXOpen::Vector3d direction, NXOpen::NXObject* body2, int boolType) {
    NXOpen::Features::Feature* nullNXOpen_Features_Feature(NULL);
    NXOpen::Features::CylinderBuilder* cylinderBuilder1;
    cylinderBuilder1 = workPart->Features()->CreateCylinderBuilder(nullNXOpen_Features_Feature);
    cylinderBuilder1->Diameter()->SetValue(diameter);
    cylinderBuilder1->Height()->SetValue(height);
    NXOpen::GeometricUtilities::BooleanOperation::BooleanType type;
    switch (boolType) {
    case 0:type = NXOpen::GeometricUtilities::BooleanOperation::BooleanTypeCreate; break;
    case 1:type = NXOpen::GeometricUtilities::BooleanOperation::BooleanTypeUnite; break;
    case 2:type = NXOpen::GeometricUtilities::BooleanOperation::BooleanTypeSubtract; break;
    case 3:type = NXOpen::GeometricUtilities::BooleanOperation::BooleanTypeIntersect; break;
    default:type = NXOpen::GeometricUtilities::BooleanOperation::BooleanTypeCreate; break;
    }
    if (boolType > 0) {
        cylinderBuilder1->BooleanOption()->SetType(type);
        NXOpen::Features::BodyFeature* bod2(dynamic_cast<NXOpen::Features::BodyFeature*>(body2));
        std::vector <NXOpen::Face*> FaceAr = bod2->GetFaces();
        std::vector<NXOpen::Body*> targetBodies2(1);
        NXOpen::Body* body1(FaceAr[0]->GetBody());
        targetBodies2[0] = body1;
        cylinderBuilder1->BooleanOption()->SetTargetBodies(targetBodies2);
    }
    NXOpen::Axis* axis1;
    axis1 = workPart->Axes()->CreateAxis(origin, direction, NXOpen::SmartObject::UpdateOptionWithinModeling);
    cylinderBuilder1->SetAxis(axis1);
    NXOpen::NXObject* cyl1;
    cyl1 = cylinderBuilder1->Commit();
    cylinderBuilder1->Destroy();
    return cyl1;
}

NXOpen::NXObject* Joint::buildCone(double Diameter, double height, double halfAngle, NXOpen::Point3d origin, NXOpen::Vector3d direction, int coneType, NXOpen::NXObject* body2, int boolType) {
    NXOpen::Features::Cone* nullNXOpen_Features_Cone(NULL);
    NXOpen::Features::ConeBuilder* coneBuilder1;
    coneBuilder1 = workPart->Features()->CreateConeBuilder(nullNXOpen_Features_Cone);
    NXOpen::Features::ConeBuilder::Types type1;
    switch (coneType) {
    case 0:type1 = NXOpen::Features::ConeBuilder::TypesTopDiameterHeightAndHalfAngle; coneBuilder1->TopDiameter()->SetValue(Diameter); break;
    case 1:type1 = NXOpen::Features::ConeBuilder::TypesBaseDiameterHeightAndHalfAngle; coneBuilder1->BaseDiameter()->SetValue(Diameter); break;
    }
    coneBuilder1->SetType(type1);
    coneBuilder1->Height()->SetValue(height);
    coneBuilder1->HalfAngle()->SetValue(halfAngle);

    NXOpen::GeometricUtilities::BooleanOperation::BooleanType type;
    switch (boolType) {
    case 0:type = NXOpen::GeometricUtilities::BooleanOperation::BooleanTypeCreate; break;
    case 1:type = NXOpen::GeometricUtilities::BooleanOperation::BooleanTypeUnite; break;
    case 2:type = NXOpen::GeometricUtilities::BooleanOperation::BooleanTypeSubtract; break;
    case 3:type = NXOpen::GeometricUtilities::BooleanOperation::BooleanTypeIntersect; break;
    default:type = NXOpen::GeometricUtilities::BooleanOperation::BooleanTypeCreate; break;
    }
    if (boolType > 0) {
        coneBuilder1->BooleanOption()->SetType(type);
        NXOpen::Features::BodyFeature* bod2(dynamic_cast<NXOpen::Features::BodyFeature*>(body2));
        std::vector <NXOpen::Face*> FaceAr = bod2->GetFaces();
        std::vector<NXOpen::Body*> targetBodies2(1);
        NXOpen::Body* body1(FaceAr[0]->GetBody());
        targetBodies2[0] = body1;
        coneBuilder1->BooleanOption()->SetTargetBodies(targetBodies2);
    }
    NXOpen::Axis* axis1;
    axis1 = workPart->Axes()->CreateAxis(origin, direction, NXOpen::SmartObject::UpdateOptionWithinModeling);
    coneBuilder1->SetAxis(axis1);

    NXOpen::NXObject* cone1;
    cone1 = coneBuilder1->Commit();
    coneBuilder1->Destroy();
    return cone1;
}

NXOpen::NXObject* Joint::buildConeByBaseDiameterHeight(double baseDiameter, double height, double halfAngle, NXOpen::Point3d origin, NXOpen::Vector3d direction, NXOpen::NXObject* body2, int boolType) {
    NXOpen::Features::Cone* nullNXOpen_Features_Cone(NULL);
    NXOpen::Features::ConeBuilder* coneBuilder1;
    coneBuilder1 = workPart->Features()->CreateConeBuilder(nullNXOpen_Features_Cone);
    coneBuilder1->SetType(NXOpen::Features::ConeBuilder::TypesBaseDiameterHeightAndHalfAngle);
    coneBuilder1->BaseDiameter()->SetValue(baseDiameter);
    coneBuilder1->Height()->SetValue(height);
    coneBuilder1->HalfAngle()->SetValue(halfAngle);
    coneBuilder1->BooleanOption()->SetType(NXOpen::GeometricUtilities::BooleanOperation::BooleanTypeCreate);

    NXOpen::Axis* axis1;
    axis1 = workPart->Axes()->CreateAxis(origin, direction, NXOpen::SmartObject::UpdateOptionWithinModeling);
    coneBuilder1->SetAxis(axis1);

    NXOpen::NXObject* cone1;
    cone1 = coneBuilder1->Commit();
    coneBuilder1->Destroy();
    return cone1;
}
NXOpen::CartesianCoordinateSystem* Joint::buildCsys(NXOpen::Point3d origin, NXOpen::Vector3d directionX, NXOpen::Vector3d directionY) {
    NXOpen::Xform* xform1;
    xform1 = workPart->Xforms()->CreateXform(origin, directionX, directionY, NXOpen::SmartObject::UpdateOptionWithinModeling, 1.0);

    NXOpen::CartesianCoordinateSystem* cartesianCoordinateSystem1;
    cartesianCoordinateSystem1 = workPart->CoordinateSystems()->CreateCoordinateSystem(xform1, NXOpen::SmartObject::UpdateOptionWithinModeling);
    return cartesianCoordinateSystem1;
}
NXOpen::NXObject* Joint::buildSketch(NXOpen::Point3d origin, NXOpen::Vector3d normal, NXOpen::Vector3d axisX) {
    NXOpen::Sketch* nullNXOpen_Sketch(NULL);
    NXOpen::SketchInPlaceBuilder* sketchInPlaceBuilder1;
    sketchInPlaceBuilder1 = workPart->Sketches()->CreateSketchInPlaceBuilder2(nullNXOpen_Sketch);

    NXOpen::Plane* plane1 = workPart->Planes()->CreatePlane(origin, normal, NXOpen::SmartObject::UpdateOption::UpdateOptionWithinModeling);
    NXOpen::Direction* dir1 = workPart->Directions()->CreateDirection(origin, axisX, NXOpen::SmartObject::UpdateOption::UpdateOptionWithinModeling);

    NXOpen::Point* po1;
    po1 = workPart->Points()->CreatePoint(origin);


    sketchInPlaceBuilder1->SetPlaneReference(plane1);
    sketchInPlaceBuilder1->SetAxisReference(dir1);
    sketchInPlaceBuilder1->SetSketchOrigin(po1);


    NXOpen::NXObject* ob1;
    ob1 = sketchInPlaceBuilder1->Commit();

    sketchInPlaceBuilder1->Destroy();

    return ob1;
}
NXOpen::NXObject* Joint::helixByLawCurve(double startAngle, double pitch, double startLimit, double endLimit, NXOpen::CartesianCoordinateSystem* csys, NXOpen::Line* baseLine, NXOpen::Line* lines[], int linesCount, bool reverse_dir) {
    NXOpen::Features::Helix* nullNXOpen_Features_Helix(NULL);
    NXOpen::Features::HelixBuilder* helixBuilder1;
    helixBuilder1 = workPart->Features()->CreateHelixBuilder(nullNXOpen_Features_Helix);
    helixBuilder1->SetOrientationOption(NXOpen::Features::HelixBuilder::OrientationOptionsSpecified);
    helixBuilder1->StartAngle()->SetValue(startAngle);
    helixBuilder1->SetSizeOption(NXOpen::Features::HelixBuilder::SizeOptionsRadius);
    helixBuilder1->SizeLaw()->SetLawType(NXOpen::GeometricUtilities::LawBuilder::TypeByLawCurve);
    helixBuilder1->PitchLaw()->Value()->SetValue(pitch);
    helixBuilder1->PitchLaw()->StartValue()->SetFormula("5");
    helixBuilder1->PitchLaw()->EndValue()->SetFormula("5");
    helixBuilder1->StartLimit()->SetPercentUsed(false);
    helixBuilder1->StartLimit()->Expression()->SetValue(startLimit);
    helixBuilder1->EndLimit()->SetPercentUsed(false);
    helixBuilder1->SetCoordinateSystem(csys);
    helixBuilder1->SizeLaw()->LawCurve()->SetAllowedEntityTypes(NXOpen::Section::AllowTypesOnlyCurves);

    for (int i = 0; i < linesCount; i++) {
        std::vector<NXOpen::IBaseCurve*> curves1(1);
        NXOpen::Line* line1 = lines[i];
        NXOpen::CurveDumbRule* curveDumbRule1;
        curveDumbRule1 = workPart->ScRuleFactory()->CreateRuleBaseCurveDumb({ lines[i] });
        helixBuilder1->SizeLaw()->LawCurve()->AllowSelfIntersection(true);

        std::vector<NXOpen::SelectionIntentRule*> rules1(1);
        rules1[0] = curveDumbRule1;
        NXOpen::NXObject* nullNXOpen_NXObject(NULL);
        NXOpen::Point3d helpPoint1(0, 0, 0);
        helixBuilder1->SizeLaw()->LawCurve()->AddToSection(rules1, lines[i], nullNXOpen_NXObject, nullNXOpen_NXObject, helpPoint1, NXOpen::Section::ModeCreate, false);
        helixBuilder1->Evaluate();
    }

    helixBuilder1->SizeLaw()->BaseLine()->SetValue(baseLine);
    /* helixBuilder1->SizeLaw()->SetReverseDirection(true);
     helixBuilder1->Evaluate();
     helixBuilder1->Evaluate();*/
    helixBuilder1->SizeLaw()->SetReverseDirection(reverse_dir);
    helixBuilder1->Evaluate();
    //helixBuilder1->EndLimit()->Expression()->SetFormula("75");
    helixBuilder1->EndLimit()->Expression()->SetValue(endLimit);
    NXOpen::NXObject* hel1;
    hel1 = helixBuilder1->Commit();
    helixBuilder1->Destroy();
    return hel1;
}

NXOpen::NXObject* Joint::subtractBodies(NXOpen::NXObject* body1, NXOpen::NXObject* body2) {

    NXOpen::Features::BooleanFeature* nullNXOpen_Features_BooleanFeature(NULL);
    NXOpen::Features::BooleanBuilder* booleanBuilder1;
    booleanBuilder1 = workPart->Features()->CreateBooleanBuilderUsingCollector(nullNXOpen_Features_BooleanFeature);

    NXOpen::ScCollector* scCollector1;
    scCollector1 = booleanBuilder1->ToolBodyCollector();

    NXOpen::GeometricUtilities::BooleanRegionSelect* booleanRegionSelect1;
    booleanRegionSelect1 = booleanBuilder1->BooleanRegionSelect();

    booleanBuilder1->SetTolerance(0.01);

    booleanBuilder1->SetOperation(NXOpen::Features::Feature::BooleanTypeSubtract);
    NXOpen::Features::BodyFeature* bod1(dynamic_cast<NXOpen::Features::BodyFeature*>(body1));
    std::vector <NXOpen::Face*> FaceAr = bod1->GetFaces();
    NXOpen::Body* sub_body1 = FaceAr[0]->GetBody();
    bool added1;
    added1 = booleanBuilder1->Targets()->Add(sub_body1);

    std::vector<NXOpen::TaggedObject*> targets1(1);
    targets1[0] = sub_body1;
    booleanRegionSelect1->AssignTargets(targets1);

    NXOpen::ScCollector* scCollector2;
    scCollector2 = workPart->ScCollectors()->CreateCollector();

    NXOpen::Features::BodyFeature* bod2(dynamic_cast<NXOpen::Features::BodyFeature*>(body2));
    std::vector <NXOpen::Face*> FaceAr2 = bod2->GetFaces();

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
    return sub1;
}



/**
* Calculate the intersection between two lines, which defined by point on line and angle
* works only with X,Y coordinates!
* */
NXOpen::Point3d Joint::getLinesIntersectionPoint(NXOpen::Point3d a, double ang1, NXOpen::Point3d b, double ang2) {
    double n1x = sin(ang1 * DEGRA);
    double n1y = -cos(ang1 * DEGRA);
    double n2x = sin(ang2 * DEGRA);
    double n2y = -cos(ang2 * DEGRA);
    double c = n1x * a.X + n1y * a.Y;
    double f = n2x * b.X + n2y * b.Y;
    NXOpen::Point3d result;
    result.Y = (c * n2x - n1x * f) / (n2x * n1y - n1x * n2y);
    result.X = (f - n2y * result.Y) / n2x;
    result.Z = 0;
    return result;
}
NXOpen::Point3d Joint::multiplyByMatrix4x4(NXOpen::Point3d p1, NXOpen::Matrix4x4 mat) {
    NXOpen::Point4d p(p1.X, p1.Y, p1.Z, 1);
    NXOpen::Point3d result;
    result.X = mat.Rxx * p.X + mat.Rxy * p.Y + mat.Rxz * p.Z + mat.Xt * p.W;
    result.Y = mat.Ryx * p.X + mat.Ryy * p.Y + mat.Ryz * p.Z + mat.Yt * p.W;
    result.Z = mat.Rzx * p.X + mat.Rzy * p.Y + mat.Rzz * p.Z + mat.Zt * p.W;
    return result;

}