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
    cone_recess_len = { 16 };//конусная выточка в муфте, длина
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


    //эскиз закона изменения диаметра спирали---------------------------------------------------------------------------------------------------------------
    NXOpen::Sketch* sk1(dynamic_cast<NXOpen::Sketch*>(buildSketch({ 0,0,0 }, { 0,0,1 }, { 1,0,0 })));
    sk1->Activate(NXOpen::Sketch::ViewReorientTrue);


    //базовая линия
    NXOpen::Point3d pl1(0, 0, 0);
    NXOpen::Point3d p2(50, 0, 0);
    NXOpen::Line* l1 = workPart->Curves()->CreateLine(pl1, p2);
    theSession->ActiveSketch()->AddGeometry(l1, NXOpen::Sketch::InferConstraintsOptionInferNoConstraints);
    //основные линии ниппеля
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



    //система координат для спирали---------------------------------------------------------------------------------------------------------------------


    NXOpen::CartesianCoordinateSystem* cartesianCoordinateSystem4 = buildCsys({ 0,0,ln }, { 1,0,0 }, { 0,-1,0 });
    //спираль по закону, заданному прямыми линиями---------------------------------------------------------------------------------------------------------------

    NXOpen::Line* hel_lines1[] = { l2,l3 };
    NXOpen::NXObject* hel1 = helixByLawCurve(90.0, P, 0.0, ln - thr_esc + 1.5 * P, cartesianCoordinateSystem4, l1, hel_lines1, 2, false);

    NXOpen::CartesianCoordinateSystem* cartesianCoordinateSystem5 = buildCsys({ 0,0,box_l1 + rotZ_Offset }, { 1,0,0 }, { 0,-1,0 });
    NXOpen::Line* hel_lines2[] = { hel2l1,hel2l2 };
    NXOpen::NXObject* hel2 = helixByLawCurve(270, P, -1.5 * P, box_l1, cartesianCoordinateSystem5, l1, hel_lines2, 2, true);


    //эскиз профиля резьбы--------------------------------------------------------------------------------------------------------




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

    //вычитание конуса и заметания

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
            face->SetColor(108);//цилинлрическая поверхность ниппеля
            count++;
        }
        if (type_face == UF_MODL_PLANAR_FACE)
        {
            UF_MODL_ask_face_parm_2(face->Tag(), point_on_face, parms, face_point);
            if (face_point[2] >= -1) face->SetColor(6);//плоская поверхность ниппеля
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
            if (face_point[2] >= 5) face->SetColor(211);//дальняя плоская поверхность муты
            else face->SetColor(83);//ближняя плоская поверхность муфты

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


    //скрытие элементов дерева построения

    std::vector<NXOpen::DisplayableObject*> dis_objects2(6);
    dis_objects2[0] = sk1;
    dis_objects2[1] = sk2;
    dis_objects2[2] = sk3;
    dis_objects2[3] = dynamic_cast<NXOpen::Spline*> (hel1);
    dis_objects2[4] = dynamic_cast<NXOpen::Spline*> (hel2);


    theSession->DisplayManager()->BlankObjects(dis_objects2);

    theSession->ApplicationSwitchImmediate("UG_APP_SFEM");



    // ----------------------------------------------
    //   Menu: Вставить->Подготовка модели->Обрезка->Разделить тело...
    // ----------------------------------------------
    NXOpen::NXObject* split_body;
    {
        NXOpen::Features::SplitBody* nullNXOpen_Features_SplitBody(NULL);
        NXOpen::Features::SplitBodyBuilder* splitBodyBuilder1;
        splitBodyBuilder1 = workPart->Features()->CreateSplitBodyBuilderUsingCollector(nullNXOpen_Features_SplitBody);

        NXOpen::Point3d origin1(0.0, 0.0, 0.0);
        NXOpen::Vector3d normal1(0.0, 0.0, 1.0);
        NXOpen::Plane* plane1;
        plane1 = workPart->Planes()->CreatePlane(origin1, normal1, NXOpen::SmartObject::UpdateOptionWithinModeling);

        splitBodyBuilder1->BooleanTool()->FacePlaneTool()->SetToolPlane(plane1);

        splitBodyBuilder1->BooleanTool()->SetToolOption(NXOpen::GeometricUtilities::BooleanToolBuilder::BooleanToolTypeNewPlane);

        splitBodyBuilder1->SetKeepImprintedEdges(true);

        splitBodyBuilder1->BooleanTool()->ExtrudeRevolveTool()->ToolSection()->PrepareMappingData();


        splitBodyBuilder1->BooleanTool()->ExtrudeRevolveTool()->ToolSection()->SetDistanceTolerance(0.01);

        splitBodyBuilder1->BooleanTool()->ExtrudeRevolveTool()->ToolSection()->SetChainingTolerance(0.0094999999999999998);

        NXOpen::ScCollector* scCollector1;
        scCollector1 = workPart->ScCollectors()->CreateCollector();

        std::vector<NXOpen::Body*> bodies1(1);
        NXOpen::Body* body1(FaceAr[0]->GetBody());
        bodies1[0] = body1;
        NXOpen::BodyDumbRule* bodyDumbRule1;
        bodyDumbRule1 = workPart->ScRuleFactory()->CreateRuleBodyDumb(bodies1, true);

        std::vector<NXOpen::SelectionIntentRule*> rules1(1);
        rules1[0] = bodyDumbRule1;
        scCollector1->ReplaceRules(rules1, false);

        splitBodyBuilder1->SetTargetBodyCollector(scCollector1);

        std::vector<NXOpen::Body*> bodies2(2);
        bodies2[0] = body1;
        NXOpen::Body* body2(FaceAr1[0]->GetBody());
        bodies2[1] = body2;
        NXOpen::BodyDumbRule* bodyDumbRule2;
        bodyDumbRule2 = workPart->ScRuleFactory()->CreateRuleBodyDumb(bodies2, true);

        std::vector<NXOpen::SelectionIntentRule*> rules2(1);
        rules2[0] = bodyDumbRule2;
        scCollector1->ReplaceRules(rules2, false);

        splitBodyBuilder1->SetTargetBodyCollector(scCollector1);

        plane1->SetMethod(NXOpen::PlaneTypes::MethodTypeDistance);

        std::vector<NXOpen::NXObject*> geom1(1);
        NXOpen::DatumPlane* datumPlane1(dynamic_cast<NXOpen::DatumPlane*>(workPart->Datums()->FindObject("DATUM_CSYS(0) YZ plane")));
        geom1[0] = datumPlane1;
        plane1->SetGeometry(geom1);

        plane1->SetFlip(false);

        plane1->SetReverseSide(false);

        plane1->SetAlternate(NXOpen::PlaneTypes::AlternateTypeOne);

        plane1->Evaluate();

        plane1->SetMethod(NXOpen::PlaneTypes::MethodTypeDistance);

        std::vector<NXOpen::NXObject*> geom2(1);
        geom2[0] = datumPlane1;
        plane1->SetGeometry(geom2);

        plane1->SetFlip(false);

        plane1->SetReverseSide(false);

        NXOpen::Expression* expression4;
        expression4 = plane1->Expression();

        expression4->SetRightHandSide("0");

        plane1->SetAlternate(NXOpen::PlaneTypes::AlternateTypeOne);

        plane1->Evaluate();

        split_body = splitBodyBuilder1->Commit();

        splitBodyBuilder1->BooleanTool()->ExtrudeRevolveTool()->ToolSection()->CleanMappingData();

        splitBodyBuilder1->Destroy();

    }



    // ----------------------------------------------
       //   Menu: Файл->Утилиты->Новый КЭ и симуляция...
       // ----------------------------------------------

    NXOpen::BasePart::Units units1;
    units1 = workPart->PartUnits();

    NXOpen::CAE::Xyplot::BaseTemplateManager* baseTemplateManager1;
    baseTemplateManager1 = theSession->XYPlotManager()->TemplateManager();

    // NXOpen::NXString path = workPart->FullPath();
    std::string path = workPart->FullPath().GetText();
    path = path.erase(path.size() - 4, 4);


    NXOpen::BasePart* basePart1;
    basePart1 = theSession->Parts()->NewBaseDisplay(path + "_fem1.fem", NXOpen::BasePart::UnitsMillimeters);

    workPart = NULL;
    NXOpen::CAE::FemPart* workFemPart(dynamic_cast<NXOpen::CAE::FemPart*>(theSession->Parts()->BaseWork()));
    displayPart = NULL;
    NXOpen::CAE::FemPart* displayFemPart(dynamic_cast<NXOpen::CAE::FemPart*>(theSession->Parts()->BaseDisplay()));
    NXOpen::CAE::FemPart* femPart1(dynamic_cast<NXOpen::CAE::FemPart*>(workFemPart));
    femPart1->PolygonGeometryMgr()->SetPolygonBodyResolutionOnFemBodies(NXOpen::CAE::PolygonGeometryManager::PolygonBodyResolutionTypeHigh);

    NXOpen::CAE::FemPart* femPart2(dynamic_cast<NXOpen::CAE::FemPart*>(workFemPart));
    NXOpen::CAE::FemCreationOptions* femCreationOptions1;
    femCreationOptions1 = femPart2->NewFemCreationOptions();

    NXOpen::CAE::FemPart* femPart3(dynamic_cast<NXOpen::CAE::FemPart*>(workFemPart));
    NXOpen::CAE::FemSynchronizeOptions* femSynchronizeOptions1;
    femSynchronizeOptions1 = femPart3->NewFemSynchronizeOptions();

    femSynchronizeOptions1->SetSynchronizePointsFlag(false);

    femSynchronizeOptions1->SetSynchronizeCreateMeshPointsFlag(false);

    femSynchronizeOptions1->SetSynchronizeCoordinateSystemFlag(false);

    femSynchronizeOptions1->SetSynchronizeLinesFlag(false);

    femSynchronizeOptions1->SetSynchronizeArcsFlag(false);

    femSynchronizeOptions1->SetSynchronizeSplinesFlag(false);

    femSynchronizeOptions1->SetSynchronizeConicsFlag(false);

    femSynchronizeOptions1->SetSynchronizeSketchCurvesFlag(false);

    femSynchronizeOptions1->SetSynchronizeDplaneFlag(false);

    // NXOpen::Part* part1(dynamic_cast<NXOpen::Part*>(theSession->Parts()->FindObject("model1")));
    // femCreationOptions1->SetCadData(part1, "D:\\doc\\SCIENCE\\auto_models\\model1.prt");

    femCreationOptions1->SetCadData(part1, path);

    std::vector<NXOpen::Body*> bodies3(4);
    NXOpen::Body* temp_body;

    NXOpen::Features::BodyFeature* split_bodies(dynamic_cast<NXOpen::Features::BodyFeature*>(split_body));
    std::vector <NXOpen::Face*> FaceAr2 = split_bodies->GetFaces();
    temp_body = FaceAr2[0]->GetBody();
    bodies3.push_back(temp_body);
    for (int i = 0; i < FaceAr2.size(); i++) {

        if (FaceAr2[i]->GetBody() != temp_body) {
            bodies3.push_back(temp_body);
        }
    }


    femCreationOptions1->SetGeometryOptions(NXOpen::CAE::FemCreationOptions::UseBodiesOptionVisibleBodies, bodies3, femSynchronizeOptions1);

    femCreationOptions1->SetSolverOptions("NX NASTRAN", "Structural", NXOpen::CAE::BaseFemPart::AxisymAbstractionTypeNone);

    std::vector<NXOpen::NXString> description1(0);
    femCreationOptions1->SetDescription(description1);

    femCreationOptions1->SetMorphingFlag(false);

    NXOpen::CoordinateSystem* nullNXOpen_CoordinateSystem(NULL);
    femCreationOptions1->SetCyclicSymmetryData(false, nullNXOpen_CoordinateSystem);

    NXOpen::CAE::FemPart* femPart4(dynamic_cast<NXOpen::CAE::FemPart*>(workFemPart));
    femPart4->FinalizeCreation(femCreationOptions1);

    delete femSynchronizeOptions1;
    delete femCreationOptions1;
    NXOpen::CAE::Xyplot::BaseTemplateManager* baseTemplateManager2;
    baseTemplateManager2 = theSession->XYPlotManager()->TemplateManager();

    NXOpen::BasePart* basePart2;
    basePart2 = theSession->Parts()->NewBaseDisplay(path + "_sim1.sim", NXOpen::BasePart::UnitsMillimeters);

    NXOpen::CAE::SimPart* workSimPart(dynamic_cast<NXOpen::CAE::SimPart*>(theSession->Parts()->BaseWork()));
    NXOpen::CAE::SimPart* displaySimPart(dynamic_cast<NXOpen::CAE::SimPart*>(theSession->Parts()->BaseDisplay()));
    NXOpen::CAE::SimPart* simPart1(dynamic_cast<NXOpen::CAE::SimPart*>(workSimPart));
    std::vector<NXOpen::NXString> description2(0);
    simPart1->FinalizeCreation(femPart4, -1, description2);

    workSimPart->ModelingViews()->WorkView()->Regenerate();




    NXOpen::CAE::SimPart* simPart2(dynamic_cast<NXOpen::CAE::SimPart*>(workSimPart));

    simSimulation1 = simPart2->Simulation();

    NXOpen::CAE::SimSolution* simSolution2;
    simSolution2 = simSimulation1->CreateSolution("NX NASTRAN", "Structural", "ADVNL 601,106", "Solution 1", NXOpen::CAE::SimSimulation::AxisymAbstractionTypeNone);

    NXOpen::CAE::PropertyTable* propertyTable2;
    propertyTable2 = simSolution2->PropertyTable();

    NXOpen::CAE::CaePart* caePart3(dynamic_cast<NXOpen::CAE::CaePart*>(workSimPart));
    NXOpen::CAE::ModelingObjectPropertyTable* modelingObjectPropertyTable3;
    modelingObjectPropertyTable3 = caePart3->ModelingObjectPropertyTables()->CreateModelingObjectPropertyTable("Bulk Data Echo Request", "NX NASTRAN - Structural", "NX NASTRAN", "Bulk Data Echo Request1", 1);

    NXOpen::CAE::CaePart* caePart4(dynamic_cast<NXOpen::CAE::CaePart*>(workSimPart));
    NXOpen::CAE::ModelingObjectPropertyTable* modelingObjectPropertyTable4;
    modelingObjectPropertyTable4 = caePart4->ModelingObjectPropertyTables()->CreateModelingObjectPropertyTable("Structural Output Requests", "NX NASTRAN - Structural", "NX NASTRAN", "Structural Output Requests1", 2);

    NXOpen::CAE::PropertyTable* structOutputTable = modelingObjectPropertyTable4->PropertyTable();
    structOutputTable->SetBooleanPropertyValue("Contact Result - Enable", true);

    simSolution2->Rename("Solution 1", false);

    NXOpen::CAE::PropertyTable* propertyTable3;
    propertyTable3 = simSolution2->PropertyTable();

    propertyTable3->SetNamedPropertyTablePropertyValue("Bulk Data Echo Request", modelingObjectPropertyTable3);

    propertyTable3->SetNamedPropertyTablePropertyValue("Output Requests", modelingObjectPropertyTable4);

    NXOpen::CAE::SimSolutionStep* simSolutionStep1(dynamic_cast<NXOpen::CAE::SimSolutionStep*>(simSolution2->FindObject("SolutionStep[Subcase - Nonlinear Implicit]")));
    simSolution2->SetActiveStep(simSolutionStep1);
    theSession->Parts()->SetWork(femPart4);

    workFemPart = dynamic_cast<NXOpen::CAE::FemPart*>(theSession->Parts()->BaseWork()); // model1_fem1

    //материал 

    NXOpen::PhysMat::PhysicalMaterialListBuilder* physicalMaterialListBuilder1;
    physicalMaterialListBuilder1 = workFemPart->MaterialManager()->PhysicalMaterials()->CreateListBlockBuilder();

    NXOpen::PhysMat::PhysicalMaterialAssignBuilder* physicalMaterialAssignBuilder1;
    physicalMaterialAssignBuilder1 = workFemPart->MaterialManager()->PhysicalMaterials()->CreateMaterialAssignBuilder();

    NXOpen::PhysicalMaterial* physicalMaterial1(dynamic_cast<NXOpen::PhysicalMaterial*>(workFemPart->MaterialManager()->PhysicalMaterials()->LoadFromNxmatmllibrary("Steel")));
    physicalMaterial1->AssignToAllBodies();

    physicalMaterialAssignBuilder1->Destroy();

    physicalMaterialListBuilder1->Destroy();

    // ----------------------------------------------
 //   Menu: Вставить->Сетка->3D тетраэдральная сетка...
 // ----------------------------------------------

    NXOpen::CAE::FEModel* fEModel1(dynamic_cast<NXOpen::CAE::FEModel*>(femPart4->BaseFEModel()));
    NXOpen::CAE::MeshManager* meshManager1(dynamic_cast<NXOpen::CAE::MeshManager*>(fEModel1->Find("MeshManager")));
    NXOpen::CAE::Mesh3d* nullNXOpen_CAE_Mesh3d(NULL);
    NXOpen::CAE::Mesh3dTetBuilder* mesh3dTetBuilder1;
    mesh3dTetBuilder1 = meshManager1->CreateMesh3dTetBuilder(nullNXOpen_CAE_Mesh3d);

    NXOpen::CAE::MeshCollector* nullNXOpen_CAE_MeshCollector(NULL);
    mesh3dTetBuilder1->ElementType()->DestinationCollector()->SetElementContainer(nullNXOpen_CAE_MeshCollector);

    mesh3dTetBuilder1->ElementType()->SetElementTypeName("CTETRA(4)");

    mesh3dTetBuilder1->ElementType()->DestinationCollector()->SetAutomaticMode(false);

    mesh3dTetBuilder1->ElementType()->DestinationCollector()->SetAutomaticMode(true);

    mesh3dTetBuilder1->PropertyTable()->SetBooleanPropertyValue("mapped mesh option bool", false);

    mesh3dTetBuilder1->PropertyTable()->SetBooleanPropertyValue("two elements through thickness bool", true);


    NXOpen::CAE::CAEBody* cAEBody1(dynamic_cast<NXOpen::CAE::CAEBody*>(workFemPart->FindObject("CAE_Body(1)")));
    //NXOpen::CAE::CAEBody* cAEBody1(dynamic_cast<NXOpen::CAE::CAEBody*>(caeboidues[0]));
    bool added1;
    added1 = mesh3dTetBuilder1->SelectionList()->Add(cAEBody1);

    NXOpen::CAE::CAEBody* cAEBody2(dynamic_cast<NXOpen::CAE::CAEBody*>(workFemPart->FindObject("CAE_Body(4)")));
    bool added2;
    added2 = mesh3dTetBuilder1->SelectionList()->Add(cAEBody2);

    NXOpen::CAE::CAEBody* cAEBody3(dynamic_cast<NXOpen::CAE::CAEBody*>(workFemPart->FindObject("CAE_Body(3)")));
    bool added3;
    added3 = mesh3dTetBuilder1->SelectionList()->Add(cAEBody3);

    NXOpen::CAE::CAEBody* cAEBody4(dynamic_cast<NXOpen::CAE::CAEBody*>(workFemPart->FindObject("CAE_Body(2)")));
    bool added4;
    added4 = mesh3dTetBuilder1->SelectionList()->Add(cAEBody4);


    mesh3dTetBuilder1->SetAutoResetOption(false);

    mesh3dTetBuilder1->ElementType()->SetElementDimension(NXOpen::CAE::ElementTypeBuilder::ElementTypeFreeSolid);

    mesh3dTetBuilder1->ElementType()->SetElementTypeName("CTETRA(4)");

    NXOpen::CAE::DestinationCollectorBuilder* destinationCollectorBuilder1;
    destinationCollectorBuilder1 = mesh3dTetBuilder1->ElementType()->DestinationCollector();

    destinationCollectorBuilder1->SetElementContainer(nullNXOpen_CAE_MeshCollector);

    destinationCollectorBuilder1->SetAutomaticMode(true);

    NXOpen::Unit* unit3(dynamic_cast<NXOpen::Unit*>(workFemPart->UnitCollection()->FindObject("MilliMeter")));
    mesh3dTetBuilder1->PropertyTable()->SetBaseScalarWithDataPropertyValue("quad mesh overall edge size", mesh_size.c_str(), unit3);

    mesh3dTetBuilder1->PropertyTable()->SetBooleanPropertyValue("mapped mesh option bool", false);

    mesh3dTetBuilder1->PropertyTable()->SetBooleanPropertyValue("multiblock cylinder option bool", false);

    mesh3dTetBuilder1->PropertyTable()->SetIntegerPropertyValue("fillet num elements", 3);

    mesh3dTetBuilder1->PropertyTable()->SetIntegerPropertyValue("num elements on cylinder circumference", 6);

    mesh3dTetBuilder1->PropertyTable()->SetBaseScalarWithDataPropertyValue("element size on cylinder height", "1", unit3);

    mesh3dTetBuilder1->PropertyTable()->SetBooleanPropertyValue("create pyramids bool", false);

    mesh3dTetBuilder1->PropertyTable()->SetIntegerPropertyValue("midnodes", 0);

    mesh3dTetBuilder1->PropertyTable()->SetBooleanPropertyValue("geometry tolerance option bool", false);

    mesh3dTetBuilder1->PropertyTable()->SetBaseScalarWithDataPropertyValue("geometry tolerance", "0", unit3);

    NXOpen::Unit* nullNXOpen_Unit(NULL);
    mesh3dTetBuilder1->PropertyTable()->SetBaseScalarWithDataPropertyValue("max jacobian", "10", nullNXOpen_Unit);

    mesh3dTetBuilder1->PropertyTable()->SetBaseScalarWithDataPropertyValue("surface mesh size variation", "50", nullNXOpen_Unit);

    mesh3dTetBuilder1->PropertyTable()->SetBaseScalarWithDataPropertyValue("volume mesh size variation", "50", nullNXOpen_Unit);

    mesh3dTetBuilder1->PropertyTable()->SetBaseScalarWithDataPropertyValue("internal mesh gradation", "1.05", nullNXOpen_Unit);

    mesh3dTetBuilder1->PropertyTable()->SetBooleanPropertyValue("internal max edge option bool", false);

    mesh3dTetBuilder1->PropertyTable()->SetBaseScalarWithDataPropertyValue("internal max edge length value", "0", unit3);

    mesh3dTetBuilder1->PropertyTable()->SetBooleanPropertyValue("two elements through thickness bool", true);

    mesh3dTetBuilder1->PropertyTable()->SetBooleanPropertyValue("mesh transition bool", false);

    mesh3dTetBuilder1->PropertyTable()->SetBooleanPropertyValue("remesh on bad quality bool", false);

    mesh3dTetBuilder1->PropertyTable()->SetBooleanPropertyValue("maximum edge length bool", false);

    mesh3dTetBuilder1->PropertyTable()->SetBaseScalarWithDataPropertyValue("maximum edge length", "1", unit3);

    mesh3dTetBuilder1->PropertyTable()->SetBaseScalarWithDataPropertyValue("small feature tolerance", "10", nullNXOpen_Unit);

    mesh3dTetBuilder1->PropertyTable()->SetBaseScalarWithDataPropertyValue("small feature value", "0.1", nullNXOpen_Unit);

    mesh3dTetBuilder1->PropertyTable()->SetIntegerPropertyValue("boundary layer element type", 3);

    mesh3dTetBuilder1->PropertyTable()->SetBooleanPropertyValue("insert blend elements", true);

    NXOpen::Unit* unit4(dynamic_cast<NXOpen::Unit*>(workFemPart->UnitCollection()->FindObject("Degrees")));
    mesh3dTetBuilder1->PropertyTable()->SetBaseScalarWithDataPropertyValue("blending angle", "90", unit4);

    mesh3dTetBuilder1->PropertyTable()->SetBaseScalarWithDataPropertyValue("sweep angle", "45", unit4);

    mesh3dTetBuilder1->PropertyTable()->SetBooleanPropertyValue("control aspect ratio", false);

    mesh3dTetBuilder1->PropertyTable()->SetBaseScalarWithDataPropertyValue("maximum exposed aspect ratio", "1000", nullNXOpen_Unit);

    mesh3dTetBuilder1->PropertyTable()->SetBooleanPropertyValue("control slender", false);

    mesh3dTetBuilder1->PropertyTable()->SetBaseScalarWithDataPropertyValue("minimum aspect ratio", "0.01", nullNXOpen_Unit);

    mesh3dTetBuilder1->PropertyTable()->SetBaseScalarWithDataPropertyValue("maximum imprint dihedral angle", "150", unit4);

    mesh3dTetBuilder1->PropertyTable()->SetBaseScalarWithDataPropertyValue("gradation rate", "10", nullNXOpen_Unit);

    mesh3dTetBuilder1->PropertyTable()->SetBaseScalarWithDataPropertyValue("smoothing distance factor", "3", nullNXOpen_Unit);

    mesh3dTetBuilder1->PropertyTable()->SetBooleanPropertyValue("all-tet boundary layer", false);

    mesh3dTetBuilder1->PropertyTable()->SetIntegerPropertyValue("dont format mesh to solver", 0);

    mesh3dTetBuilder1->PropertyTable()->SetBaseScalarWithDataPropertyValue("quad mesh edge match tolerance", "0.02", nullNXOpen_Unit);

    mesh3dTetBuilder1->PropertyTable()->SetBaseScalarWithDataPropertyValue("quad mesh smoothness tolerance", "0.01", unit3);

    mesh3dTetBuilder1->PropertyTable()->SetBaseScalarWithDataPropertyValue("min face angle", "20", unit4);

    mesh3dTetBuilder1->PropertyTable()->SetIntegerPropertyValue("mesh time stamp", 0);

    mesh3dTetBuilder1->PropertyTable()->SetBaseScalarWithDataPropertyValue("quad mesh node coincidence tolerance", "0.0001", unit3);

    mesh3dTetBuilder1->PropertyTable()->SetIntegerPropertyValue("mesh edit allowed", 0);

    mesh3dTetBuilder1->PropertyTable()->SetBaseScalarWithDataPropertyValue("edge angle", "15", unit4);

    mesh3dTetBuilder1->PropertyTable()->SetIntegerPropertyValue("merge edge toggle", 0);

    mesh3dTetBuilder1->PropertyTable()->SetIntegerPropertyValue("auto constraining", 1);

    mesh3dTetBuilder1->PropertyTable()->SetIntegerPropertyValue("curvature scaling", 1);

    mesh3dTetBuilder1->PropertyTable()->SetBaseScalarWithDataPropertyValue("target angle", "45", unit4);

    mesh3dTetBuilder1->PropertyTable()->SetIntegerPropertyValue("edge shape", 2);

    std::vector<NXOpen::CAE::Mesh*> meshes1;
    meshes1 = mesh3dTetBuilder1->CommitMesh();


    mesh3dTetBuilder1->Destroy();


    // ----------------------------------------------
//   Menu: Вставить->Сетка->1D соединения...
// ----------------------------------------------


    std::vector<NXOpen::DisplayableObject*> mesh1d_objs;
    std::vector<NXOpen::DisplayableObject*> contact_surface1_objs;
    std::vector<NXOpen::DisplayableObject*> contact_surface2_objs;
    std::vector<NXOpen::DisplayableObject*> constraint_objs;
    //NXOpen::NXObject* mesh1d;
    //{
    //    NXOpen::CAE::FEModel* fEModel1(dynamic_cast<NXOpen::CAE::FEModel*>(workFemPart->FindObject("FEModel")));
    //    NXOpen::CAE::CAEConnection* nullNXOpen_CAE_CAEConnection(NULL);
    //    NXOpen::CAE::CAEConnectionBuilder* cAEConnectionBuilder1;
    //    cAEConnectionBuilder1 = fEModel1->CaeConnections()->CreateConnectionBuilder(nullNXOpen_CAE_CAEConnection);

    //    NXOpen::CAE::MeshCollector* nullNXOpen_CAE_MeshCollector(NULL);
    //    cAEConnectionBuilder1->ElementType()->DestinationCollector()->SetElementContainer(nullNXOpen_CAE_MeshCollector);

    //    cAEConnectionBuilder1->ElementTypeRbe3()->DestinationCollector()->SetElementContainer(nullNXOpen_CAE_MeshCollector);

    //    cAEConnectionBuilder1->ElementType()->SetElementDimension(NXOpen::CAE::ElementTypeBuilder::ElementTypeConnection);

    //    cAEConnectionBuilder1->ElementTypeRbe3()->SetElementDimension(NXOpen::CAE::ElementTypeBuilder::ElementTypeSpider);

    //    NXOpen::CAE::MeshManager* meshManager1(dynamic_cast<NXOpen::CAE::MeshManager*>(fEModel1->Find("MeshManager")));
    //    NXOpen::CAE::MeshCollectorBuilder* meshCollectorBuilder1;
    //    meshCollectorBuilder1 = meshManager1->CreateCollectorBuilder(nullNXOpen_CAE_MeshCollector, "Rigid Link Collector");

    //    meshCollectorBuilder1->SetCollectorName("RBE2 Collector(1)");

    //    NXOpen::NXObject* nXObject2;
    //    nXObject2 = meshCollectorBuilder1->Commit();


    //    meshCollectorBuilder1->Destroy();

    //    NXOpen::CAE::MeshCollector* meshCollector1(dynamic_cast<NXOpen::CAE::MeshCollector*>(nXObject2));
    //    cAEConnectionBuilder1->ElementType()->DestinationCollector()->SetElementContainer(meshCollector1);


    //    cAEConnectionBuilder1->ElementTypeRbe3()->DestinationCollector()->SetElementContainer(meshCollector1);


    //    cAEConnectionBuilder1->SetMidNode(true);

    //    cAEConnectionBuilder1->SetType(NXOpen::CAE::CAEConnectionBuilder::ConnectionTypeEnumPointToFace);

    //    cAEConnectionBuilder1->ElementType()->SetElementTypeName("RBE2");

    //    cAEConnectionBuilder1->ElementType()->DestinationCollector()->SetElementContainer(meshCollector1);

    //    cAEConnectionBuilder1->ElementTypeRbe3()->SetElementTypeName("RBE2");

    //    cAEConnectionBuilder1->SetLabel(2084855);

    //    cAEConnectionBuilder1->ElementType()->SetElementTypeName("RBE2");


    //    NXOpen::Unit* unit1(dynamic_cast<NXOpen::Unit*>(workFemPart->UnitCollection()->FindObject("MilliMeter")));
    //    NXOpen::Expression* expression1;
    //    expression1 = workFemPart->Expressions()->CreateSystemExpressionWithUnits("0", unit1);

    //    NXOpen::Expression* expression2;
    //    expression2 = workFemPart->Expressions()->CreateSystemExpressionWithUnits("p1_x=0.00000000000", unit1);

    //    NXOpen::Expression* expression3;
    //    expression3 = workFemPart->Expressions()->CreateSystemExpressionWithUnits("p2_y=0.00000000000", unit1);

    //    NXOpen::Expression* expression4;
    //    expression4 = workFemPart->Expressions()->CreateSystemExpressionWithUnits("p3_z=0.00000000000", unit1);

    //    NXOpen::Expression* expression5;
    //    expression5 = workFemPart->Expressions()->CreateSystemExpressionWithUnits("p4_xdelta=0.00000000000", unit1);

    //    NXOpen::Expression* expression6;
    //    expression6 = workFemPart->Expressions()->CreateSystemExpressionWithUnits("p5_ydelta=0.00000000000", unit1);

    //    NXOpen::Expression* expression7;
    //    expression7 = workFemPart->Expressions()->CreateSystemExpressionWithUnits("p6_zdelta=0.00000000000", unit1);

    //    NXOpen::Expression* expression8;
    //    expression8 = workFemPart->Expressions()->CreateSystemExpressionWithUnits("p7_radius=0.00000000000", unit1);

    //    NXOpen::Unit* unit2(dynamic_cast<NXOpen::Unit*>(workFemPart->UnitCollection()->FindObject("Degrees")));
    //    NXOpen::Expression* expression9;
    //    expression9 = workFemPart->Expressions()->CreateSystemExpressionWithUnits("p8_angle=0.00000000000", unit2);

    //    NXOpen::Expression* expression10;
    //    expression10 = workFemPart->Expressions()->CreateSystemExpressionWithUnits("p9_zdelta=0.00000000000", unit1);

    //    NXOpen::Expression* expression11;
    //    expression11 = workFemPart->Expressions()->CreateSystemExpressionWithUnits("p10_radius=0.00000000000", unit1);

    //    NXOpen::Expression* expression12;
    //    expression12 = workFemPart->Expressions()->CreateSystemExpressionWithUnits("p11_angle1=0.00000000000", unit2);

    //    NXOpen::Expression* expression13;
    //    expression13 = workFemPart->Expressions()->CreateSystemExpressionWithUnits("p12_angle2=0.00000000000", unit2);

    //    NXOpen::Expression* expression14;
    //    expression14 = workFemPart->Expressions()->CreateSystemExpressionWithUnits("p13_distance=0", unit1);

    //    NXOpen::Expression* expression15;
    //    expression15 = workFemPart->Expressions()->CreateSystemExpressionWithUnits("p14_arclen=0", unit1);

    //    NXOpen::Unit* nullNXOpen_Unit(NULL);
    //    NXOpen::Expression* expression16;
    //    expression16 = workFemPart->Expressions()->CreateSystemExpressionWithUnits("p15_percent=0", nullNXOpen_Unit);

    //    expression2->SetFormula("0");

    //    expression3->SetFormula("0");

    //    expression4->SetFormula("-50");

    //    expression5->SetFormula("0");

    //    expression6->SetFormula("0");

    //    expression7->SetFormula("0");

    //    expression8->SetFormula("0");

    //    expression9->SetFormula("0");

    //    expression10->SetFormula("0");

    //    expression11->SetFormula("0");

    //    expression12->SetFormula("0");

    //    expression13->SetFormula("0");

    //    expression14->SetFormula("0");

    //    expression16->SetFormula("100");

    //    expression15->SetFormula("0");


    //    NXOpen::Expression* expression17;
    //    expression17 = workFemPart->Expressions()->CreateSystemExpressionWithUnits("p16_x=0.00000000000", unit1);

    //    NXOpen::Scalar* scalar1;
    //    scalar1 = workFemPart->Scalars()->CreateScalarExpression(expression17, NXOpen::Scalar::DimensionalityTypeNone, NXOpen::SmartObject::UpdateOptionAfterModeling);

    //    NXOpen::Expression* expression18;
    //    expression18 = workFemPart->Expressions()->CreateSystemExpressionWithUnits("p17_y=0.00000000000", unit1);

    //    NXOpen::Scalar* scalar2;
    //    scalar2 = workFemPart->Scalars()->CreateScalarExpression(expression18, NXOpen::Scalar::DimensionalityTypeNone, NXOpen::SmartObject::UpdateOptionAfterModeling);

    //    NXOpen::Expression* expression19;
    //    expression19 = workFemPart->Expressions()->CreateSystemExpressionWithUnits("p18_z=0.00000000000", unit1);

    //    NXOpen::Scalar* scalar3;
    //    scalar3 = workFemPart->Scalars()->CreateScalarExpression(expression19, NXOpen::Scalar::DimensionalityTypeNone, NXOpen::SmartObject::UpdateOptionAfterModeling);

    //    NXOpen::Point* point1;
    //    point1 = workFemPart->Points()->CreatePoint(scalar1, scalar2, scalar3, NXOpen::SmartObject::UpdateOptionAfterModeling);

    //    expression4->SetFormula("0");

    //    expression2->SetFormula("0.00000000000");

    //    expression3->SetFormula("0.00000000000");

    //    expression4->SetFormula("0.00000000000");

    //    expression2->SetFormula("0");

    //    expression3->SetFormula("0");

    //    expression4->SetFormula("0");

    //    expression2->SetFormula("0.00000000000");

    //    expression3->SetFormula("0.00000000000");

    //    expression4->SetFormula("0.00000000000");

    //    expression5->SetFormula("0.00000000000");

    //    expression6->SetFormula("0.00000000000");

    //    expression7->SetFormula("0.00000000000");

    //    expression8->SetFormula("0.00000000000");

    //    expression9->SetFormula("0.00000000000");

    //    expression10->SetFormula("0.00000000000");

    //    expression11->SetFormula("0.00000000000");

    //    expression12->SetFormula("0.00000000000");

    //    expression13->SetFormula("0.00000000000");

    //    expression16->SetFormula("100.00000000000");

    //    // ----------------------------------------------
    //    //   Dialog Begin Point
    //    // ----------------------------------------------
    //    expression4->SetFormula("-50");

    //    workFemPart->Points()->DeletePoint(point1);

    //    expression2->SetRightHandSide("0.00000000000");

    //    expression3->SetRightHandSide("0.00000000000");

    //    expression4->SetRightHandSide("-50");

    //    NXOpen::Expression* expression20;
    //    expression20 = workFemPart->Expressions()->CreateSystemExpressionWithUnits("p3_x=0.00000000000", unit1);

    //    NXOpen::Scalar* scalar4;
    //    scalar4 = workFemPart->Scalars()->CreateScalarExpression(expression20, NXOpen::Scalar::DimensionalityTypeNone, NXOpen::SmartObject::UpdateOptionAfterModeling);

    //    NXOpen::Expression* expression21;
    //    expression21 = workFemPart->Expressions()->CreateSystemExpressionWithUnits("p4_y=0.00000000000", unit1);

    //    NXOpen::Scalar* scalar5;
    //    scalar5 = workFemPart->Scalars()->CreateScalarExpression(expression21, NXOpen::Scalar::DimensionalityTypeNone, NXOpen::SmartObject::UpdateOptionAfterModeling);

    //    NXOpen::Expression* expression22;
    //    expression22 = workFemPart->Expressions()->CreateSystemExpressionWithUnits("p5_z=-50", unit1);

    //    NXOpen::Scalar* scalar6;
    //    scalar6 = workFemPart->Scalars()->CreateScalarExpression(expression22, NXOpen::Scalar::DimensionalityTypeNone, NXOpen::SmartObject::UpdateOptionAfterModeling);

    //    NXOpen::Point* point2;
    //    point2 = workFemPart->Points()->CreatePoint(scalar4, scalar5, scalar6, NXOpen::SmartObject::UpdateOptionAfterModeling);



    //    expression2->SetRightHandSide("0.00000000000");

    //    expression3->SetRightHandSide("0.00000000000");

    //    expression4->SetRightHandSide("-50");

    //    workFemPart->Points()->DeletePoint(point2);

    //    NXOpen::Expression* expression23;
    //    expression23 = workFemPart->Expressions()->CreateSystemExpressionWithUnits("p3_x=0.00000000000", unit1);

    //    NXOpen::Scalar* scalar7;
    //    scalar7 = workFemPart->Scalars()->CreateScalarExpression(expression23, NXOpen::Scalar::DimensionalityTypeNone, NXOpen::SmartObject::UpdateOptionAfterModeling);

    //    NXOpen::Expression* expression24;
    //    expression24 = workFemPart->Expressions()->CreateSystemExpressionWithUnits("p4_y=0.00000000000", unit1);

    //    NXOpen::Scalar* scalar8;
    //    scalar8 = workFemPart->Scalars()->CreateScalarExpression(expression24, NXOpen::Scalar::DimensionalityTypeNone, NXOpen::SmartObject::UpdateOptionAfterModeling);

    //    NXOpen::Expression* expression25;
    //    expression25 = workFemPart->Expressions()->CreateSystemExpressionWithUnits("p5_z=-50", unit1);

    //    NXOpen::Scalar* scalar9;
    //    scalar9 = workFemPart->Scalars()->CreateScalarExpression(expression25, NXOpen::Scalar::DimensionalityTypeNone, NXOpen::SmartObject::UpdateOptionAfterModeling);

    //    NXOpen::Point* point3;
    //    point3 = workFemPart->Points()->CreatePoint(scalar7, scalar8, scalar9, NXOpen::SmartObject::UpdateOptionAfterModeling);


    //    try
    //    {
    //        // Выражение используется
    //        workFemPart->Expressions()->Delete(expression2);
    //    }
    //    catch (const NXOpen::NXException& ex)
    //    {
    //        ex.AssertErrorCode(1050029);
    //    }

    //    try
    //    {
    //        // Выражение используется
    //        workFemPart->Expressions()->Delete(expression3);
    //    }
    //    catch (const NXOpen::NXException& ex)
    //    {
    //        ex.AssertErrorCode(1050029);
    //    }

    //    try
    //    {
    //        // Выражение используется
    //        workFemPart->Expressions()->Delete(expression4);
    //    }
    //    catch (const NXOpen::NXException& ex)
    //    {
    //        ex.AssertErrorCode(1050029);
    //    }

    //    try
    //    {
    //        // Выражение используется
    //        workFemPart->Expressions()->Delete(expression5);
    //    }
    //    catch (const NXOpen::NXException& ex)
    //    {
    //        ex.AssertErrorCode(1050029);
    //    }

    //    try
    //    {
    //        // Выражение используется
    //        workFemPart->Expressions()->Delete(expression6);
    //    }
    //    catch (const NXOpen::NXException& ex)
    //    {
    //        ex.AssertErrorCode(1050029);
    //    }

    //    try
    //    {
    //        // Выражение используется
    //        workFemPart->Expressions()->Delete(expression7);
    //    }
    //    catch (const NXOpen::NXException& ex)
    //    {
    //        ex.AssertErrorCode(1050029);
    //    }

    //    try
    //    {
    //        // Выражение используется
    //        workFemPart->Expressions()->Delete(expression8);
    //    }
    //    catch (const NXOpen::NXException& ex)
    //    {
    //        ex.AssertErrorCode(1050029);
    //    }

    //    try
    //    {
    //        // Выражение используется
    //        workFemPart->Expressions()->Delete(expression9);
    //    }
    //    catch (const NXOpen::NXException& ex)
    //    {
    //        ex.AssertErrorCode(1050029);
    //    }

    //    try
    //    {
    //        // Выражение используется
    //        workFemPart->Expressions()->Delete(expression10);
    //    }
    //    catch (const NXOpen::NXException& ex)
    //    {
    //        ex.AssertErrorCode(1050029);
    //    }

    //    try
    //    {
    //        // Выражение используется
    //        workFemPart->Expressions()->Delete(expression11);
    //    }
    //    catch (const NXOpen::NXException& ex)
    //    {
    //        ex.AssertErrorCode(1050029);
    //    }

    //    try
    //    {
    //        // Выражение используется
    //        workFemPart->Expressions()->Delete(expression12);
    //    }
    //    catch (const NXOpen::NXException& ex)
    //    {
    //        ex.AssertErrorCode(1050029);
    //    }

    //    try
    //    {
    //        // Выражение используется
    //        workFemPart->Expressions()->Delete(expression13);
    //    }
    //    catch (const NXOpen::NXException& ex)
    //    {
    //        ex.AssertErrorCode(1050029);
    //    }

    //    try
    //    {
    //        // Выражение используется
    //        workFemPart->Expressions()->Delete(expression14);
    //    }
    //    catch (const NXOpen::NXException& ex)
    //    {
    //        ex.AssertErrorCode(1050029);
    //    }

    //    try
    //    {
    //        // Выражение используется
    //        workFemPart->Expressions()->Delete(expression15);
    //    }
    //    catch (const NXOpen::NXException& ex)
    //    {
    //        ex.AssertErrorCode(1050029);
    //    }

    //    try
    //    {
    //        // Выражение используется
    //        workFemPart->Expressions()->Delete(expression16);
    //    }
    //    catch (const NXOpen::NXException& ex)
    //    {
    //        ex.AssertErrorCode(1050029);
    //    }

    //    workFemPart->Expressions()->Delete(expression1);


    //    NXOpen::Point3d point4(0.0, 0.0, -50.0);
    //    bool added1;
    //    added1 = cAEConnectionBuilder1->SourceSelection()->Add(point3, displaySimPart->ModelingViews()->WorkView(), point4);


    //    //синхронизация CAD-свойств
    //    NXOpen::CAE::BodyCollection* caeBodies(workFemPart->Bodies());

    //    for (auto it = caeBodies->begin(); it != caeBodies->end(); it++)
    //    {
    //        (*it)->SynchronizeCadProperties();
    //    }

    //    //выбор граней по цвету
    //    tag_t objTag = NULL_TAG;
    //    NXOpen::DisplayableObject* obj;
    //    int obj_type;
    //    int obj_subtype;


    //    UF_OBJ_cycle_objs_in_part(workFemPart->Tag(), UF_caegeom_type, &objTag);

    //    do {
    //        UF_OBJ_ask_type_and_subtype(objTag, &obj_type, &obj_subtype);
    //        if (obj_subtype == UF_caegeom_face_subtype) {
    //            obj = dynamic_cast<NXOpen::DisplayableObject*> (NXOpen::NXObjectManager::Get(objTag));
    //            if (obj->Color() == 108) {
    //                mesh1d_objs.push_back(obj);
    //            }




    //            //UF_MODL_ask_face_type(objTag, &type_face);
    //            //if (type_face == UF_MODL_CYLINDRICAL_FACE) i++;
    //        }
    //        UF_OBJ_cycle_objs_in_part(workFemPart->Tag(), UF_caegeom_type, &objTag);
    //    } while (objTag != NULL_TAG);


    //    bool added2;
    //    added2 = cAEConnectionBuilder1->TargetSelection()->Add(mesh1d_objs[0]);


    //    bool added3;
    //    added3 = cAEConnectionBuilder1->TargetSelection()->Add(mesh1d_objs[1]);


    //    cAEConnectionBuilder1->ElementType()->SetElementTypeName("RBE2");

    //    cAEConnectionBuilder1->SetNodeFaceProximity(0.0);

    //    cAEConnectionBuilder1->SetSearchDistance(82.462113193446598);


    //    mesh1d = cAEConnectionBuilder1->Commit();


    //    cAEConnectionBuilder1->Destroy();     }

    //simulation

              // NXOpen::CAE::SimPart* simPart1(dynamic_cast<NXOpen::CAE::SimPart*>(displaySimPart));
    NXOpen::PartLoadStatus* partLoadStatus1;
    NXOpen::PartCollection::SdpsStatus status1;
    status1 = theSession->Parts()->SetActiveDisplay(simPart1, NXOpen::DisplayPartOptionAllowAdditional, NXOpen::PartDisplayPartWorkPartOptionSameAsDisplay, &partLoadStatus1);

    // NXOpen::CAE::SimPart* workSimPart(dynamic_cast<NXOpen::CAE::SimPart*>(theSession->Parts()->BaseWork()));
    delete partLoadStatus1;
    NXOpen::CAE::CaePart* caePart1(dynamic_cast<NXOpen::CAE::CaePart*>(workSimPart));

    //NXOpen::CAE::SimSimulation* simSimulation1;
    simSimulation1 = simPart1->Simulation();

    NXOpen::CAE::SimBCBuilder* simBCBuilder1;
    simBCBuilder1 = simSimulation1->CreateBcBuilderForSimulationObjectDescriptor("Manual Surface to Surface Contact", "Face Contact(1)", 1);

    NXOpen::CAE::PropertyTable* contact_property;
    contact_property = simBCBuilder1->PropertyTable();

    NXOpen::CAE::SetManager* setManager1;
    setManager1 = simBCBuilder1->TargetSetManager();



    NXOpen::CAE::CaeRegion* nullNXOpen_CAE_CaeRegion(NULL);
    NXOpen::CAE::CaeRegionBuilder* caeRegionBuilder1;
    caeRegionBuilder1 = simSimulation1->CreateCaeRegionBuilder("Region", nullNXOpen_CAE_CaeRegion);

    NXOpen::CAE::SetManager* setManagerRegion1;
    setManagerRegion1 = caeRegionBuilder1->TargetSetManager();

    tag_t objTag = NULL_TAG;
    NXOpen::DisplayableObject* obj;
    int obj_type;
    int obj_subtype;
    std::vector<NXOpen::DisplayableObject*> postResultSurfaces;
    UF_OBJ_cycle_objs_in_part(workSimPart->Tag(), UF_caegeom_type, &objTag);

    do {
        UF_OBJ_ask_type_and_subtype(objTag, &obj_type, &obj_subtype);
        if (obj_subtype == UF_caegeom_face_subtype) {
            obj = dynamic_cast<NXOpen::DisplayableObject*> (NXOpen::NXObjectManager::Get(objTag));

            if ((obj->Color() == 6) or (obj->Color() == 36)) {
                if (obj->Color() == 6) postResultSurfaces.push_back(obj);
                contact_surface1_objs.push_back(obj);
            }
            if (obj->Color() == 83) {
                contact_surface2_objs.push_back(obj);
            }
            if (obj->Color() == 211) {
                constraint_objs.push_back(obj);
            }



            //UF_MODL_ask_face_type(objTag, &type_face);
            //if (type_face == UF_MODL_CYLINDRICAL_FACE) i++;
        }
        UF_OBJ_cycle_objs_in_part(workSimPart->Tag(), UF_caegeom_type, &objTag);
    } while (objTag != NULL_TAG);


    NXOpen::UI* theUI = NXOpen::UI::GetUI();

    //выделение узлов КЭ сетки по граням
    std::vector <NXOpen::CAE::FENode*> nodesObject1;
    std::vector <NXOpen::CAE::FENode*> nodesObject2;

    // std::vector <int> intnodesObject2;
      //  NXOpen::CAE::FEElemFace* postResultSurfacesElement(dynamic_cast<NXOpen::CAE::FEElemFace*>(postResultSurfaces[0]));
       // postResultSurfacesElement->GetElementsAndFaceIds(intnodesObject2);
       // theUI->NXMessageBox()->Show("face_ids", NXOpen::NXMessageBox::DialogTypeInformation, std::to_string(intnodesObject2[0]));
        // nodesObject1 = postResultSurfacesElement->GetNodes();

       //  NXOpen::CAE::FEElement* postResultSurfacesElement2(dynamic_cast<NXOpen::CAE::FEElement*>(postResultSurfaces[1]));
        // nodesObject2 = postResultSurfacesElement2->GetNodes();
        // NXOpen::CAE::FEElemFace* fdf;
        // fdf->GetElementsAndFaceIds();
    std::vector<NXOpen::CAE::CAEFace*> seeds1(2);
    seeds1[0] = dynamic_cast<NXOpen::CAE::CAEFace*>(postResultSurfaces[0]);
    seeds1[1] = dynamic_cast<NXOpen::CAE::CAEFace*>(postResultSurfaces[1]);
    NXOpen::CAE::RelatedNodeMethod* relatedNodeMethod1;
    relatedNodeMethod1 = caePart1->SmartSelectionMgr()->CreateRelatedNodeMethod(seeds1, true);
    nodesObject1 = relatedNodeMethod1->GetNodes();
    /*
    string nodStr;
    for (int i = 0; i < 3; i++) {
        nodStr += std::to_string(nodesObject1[i]->Label());
    }*/
    // theUI->NXMessageBox()->Show("face_ids", NXOpen::NXMessageBox::DialogTypeInformation, std::to_string(nodesObject1.size()));

    int counter = 0;
    for (auto& face : contact_surface1_objs) {

        counter++;
    }
    std::vector<NXOpen::CAE::SetObject> objects1(counter);
    for (int i = 0; i < counter; i++) {

        NXOpen::CAE::CAEFace* region_face(dynamic_cast<NXOpen::CAE::CAEFace*>(contact_surface1_objs[i]));
        //objects1.push_back({ region_face, NXOpen::CAE::CaeSetObjectSubTypeNone,counter });
        objects1[i].Obj = region_face;
        objects1[i].SubType = NXOpen::CAE::CaeSetObjectSubTypeNone;
        objects1[i].SubId = 0;

    }






    setManagerRegion1->SetTargetSetMembers(0, NXOpen::CAE::CaeSetGroupFilterTypeGeomFace, objects1);

    NXOpen::NXObject* region1;
    region1 = caeRegionBuilder1->Commit();




    NXOpen::CAE::CaeRegionBuilder* caeRegionBuilder2;
    caeRegionBuilder2 = simSimulation1->CreateCaeRegionBuilder("Region", nullNXOpen_CAE_CaeRegion);
    std::vector<NXOpen::CAE::SetObject> objects2;

    NXOpen::CAE::SetManager* setManagerRegion2;
    setManagerRegion2 = caeRegionBuilder2->TargetSetManager();

    counter = 0;
    for (auto& face : contact_surface2_objs) {

        NXOpen::CAE::CAEFace* cAEFace1(dynamic_cast<NXOpen::CAE::CAEFace*>(face));
        objects2.push_back({ cAEFace1, NXOpen::CAE::CaeSetObjectSubTypeNone,0 });

    }
    setManagerRegion2->SetTargetSetMembers(0, NXOpen::CAE::CaeSetGroupFilterTypeGeomFace, objects2);

    NXOpen::NXObject* region2;
    region2 = caeRegionBuilder2->Commit();

    NXOpen::CAE::ModelingObjectPropertyTable* modelingObjectPropertyTable1;
    modelingObjectPropertyTable1 = caePart1->ModelingObjectPropertyTables()->CreateModelingObjectPropertyTable("Contact Set Parameters", "NX NASTRAN - Structural", "NX NASTRAN", "Contact Parameters - Advanced Nonlinear Pair1", 3);

    NXOpen::CAE::PropertyTable* propertyTableAdvNonLin;
    propertyTableAdvNonLin = modelingObjectPropertyTable1->PropertyTable();

    NXOpen::Unit* second(dynamic_cast<NXOpen::Unit*>(workSimPart->UnitCollection()->FindObject("Second")));
    propertyTableAdvNonLin->SetBaseScalarWithDataPropertyValue("TZPENE", "1", second);

    NXOpen::Unit* cfactor_unit(dynamic_cast<NXOpen::Unit*>(workSimPart->UnitCollection()->FindObject("MilliMeterCubedPerMilliNewton")));
    propertyTableAdvNonLin->SetBaseScalarWithDataPropertyValue("CFACTOR1", cfactor.str(), cfactor_unit);
    propertyTableAdvNonLin->SetIntegerPropertyValue("SEGNORM", -1);





    NXOpen::CAE::CaeRegion* caeRegion1(dynamic_cast<NXOpen::CAE::CaeRegion*>(region1));
    contact_property->SetReferencePropertyValue("Source Region", caeRegion1);

    NXOpen::CAE::CaeRegion* caeRegion2(dynamic_cast<NXOpen::CAE::CaeRegion*>(region2));
    contact_property->SetReferencePropertyValue("Target Region", caeRegion2);

    NXOpen::Fields::ScalarFieldWrapper* scalarFieldWrapper1;
    scalarFieldWrapper1 = contact_property->GetScalarFieldWrapperPropertyValue("static friction coefficient");

    NXOpen::Expression* expression2;
    expression2 = scalarFieldWrapper1->GetExpression();

    expression2->SetRightHandSide("0.13");

    scalarFieldWrapper1->SetExpression(expression2);

    contact_property->SetScalarFieldWrapperPropertyValue("static friction coefficient", scalarFieldWrapper1);

    NXOpen::Fields::Field* field1;
    field1 = scalarFieldWrapper1->GetField();

    NXOpen::Fields::ScalarFieldWrapper* scalarFieldWrapper2;
    scalarFieldWrapper2 = contact_property->GetScalarFieldWrapperPropertyValue("Min Search Distance");

    NXOpen::Expression* expression3;
    expression3 = scalarFieldWrapper2->GetExpression();

    NXOpen::Unit* millimeter(dynamic_cast<NXOpen::Unit*>(workSimPart->UnitCollection()->FindObject("MilliMeter")));
    workSimPart->Expressions()->EditWithUnits(expression3, millimeter, "-0.3");

    scalarFieldWrapper2->SetExpression(expression3);

    contact_property->SetScalarFieldWrapperPropertyValue("Min Search Distance", scalarFieldWrapper2);

    NXOpen::Fields::Field* field2;
    field2 = scalarFieldWrapper2->GetField();

    NXOpen::Fields::ScalarFieldWrapper* scalarFieldWrapper3;
    scalarFieldWrapper3 = contact_property->GetScalarFieldWrapperPropertyValue("Max Search Distance");

    NXOpen::Expression* expression4;
    expression4 = scalarFieldWrapper3->GetExpression();

    workSimPart->Expressions()->EditWithUnits(expression4, millimeter, "0.3");

    scalarFieldWrapper3->SetExpression(expression4);

    contact_property->SetScalarFieldWrapperPropertyValue("Max Search Distance", scalarFieldWrapper3);

    NXOpen::Fields::Field* field3;
    field3 = scalarFieldWrapper3->GetField();

    contact_property->SetNamedPropertyTablePropertyValue("Contact Set Parameters", modelingObjectPropertyTable1);

    std::vector<NXOpen::NXString> propertyValue1(0);
    contact_property->SetTextPropertyValue("description", propertyValue1);

    NXOpen::CAE::SimLbcFolder* nullNXOpen_CAE_SimLbcFolder(NULL);
    simBCBuilder1->SetDestinationFolder(nullNXOpen_CAE_SimLbcFolder);

    NXOpen::CAE::SimBC* simBC1;
    simBC1 = simBCBuilder1->CommitAddBc();

    simBCBuilder1->Destroy();

    //заделка
    NXOpen::CAE::SimBCBuilder* fixedBuilder;
    fixedBuilder = simSimulation1->CreateBcBuilderForConstraintDescriptor("fixedConstraint", "Fixed(1)", 1);

    NXOpen::CAE::PropertyTable* propertyTable4;
    propertyTable4 = fixedBuilder->PropertyTable();

    NXOpen::CAE::SetManager* setManagerFixed;
    setManagerFixed = fixedBuilder->TargetSetManager();

    NXOpen::Fields::FieldExpression* fieldExpression1;
    fieldExpression1 = propertyTable4->GetScalarFieldPropertyValue("DOF1");

    NXOpen::Fields::FieldExpression* fieldExpression2;
    fieldExpression2 = propertyTable4->GetScalarFieldPropertyValue("DOF2");

    NXOpen::Fields::FieldExpression* fieldExpression3;
    fieldExpression3 = propertyTable4->GetScalarFieldPropertyValue("DOF3");

    NXOpen::Fields::FieldExpression* fieldExpression4;
    fieldExpression4 = propertyTable4->GetScalarFieldPropertyValue("DOF4");

    NXOpen::Fields::FieldExpression* fieldExpression5;
    fieldExpression5 = propertyTable4->GetScalarFieldPropertyValue("DOF5");

    NXOpen::Fields::FieldExpression* fieldExpression6;
    fieldExpression6 = propertyTable4->GetScalarFieldPropertyValue("DOF6");


    // ----------------------------------------------
    //   Dialog Begin Fixed Constraint
    // ----------------------------------------------
    NXOpen::Session::UndoMarkId markId31;
    markId31 = theSession->SetUndoMark(NXOpen::Session::MarkVisibilityInvisible, NXOpen::NXString("\320\227\320\260\320\264\320\265\320\273\320\272\320\260", NXOpen::NXString::UTF8));

    theSession->DeleteUndoMark(markId31, NULL);

    NXOpen::Session::UndoMarkId markId32;
    markId32 = theSession->SetUndoMark(NXOpen::Session::MarkVisibilityInvisible, NXOpen::NXString("\320\227\320\260\320\264\320\265\320\273\320\272\320\260", NXOpen::NXString::UTF8));

    std::vector<NXOpen::CAE::SetObject> objects9(2);
    NXOpen::CAE::CAEFace* cAEFace101(dynamic_cast<NXOpen::CAE::CAEFace*>(constraint_objs[0]));
    objects9[0].Obj = cAEFace101;
    objects9[0].SubType = NXOpen::CAE::CaeSetObjectSubTypeNone;
    objects9[0].SubId = 0;
    NXOpen::CAE::CAEFace* cAEFace102(dynamic_cast<NXOpen::CAE::CAEFace*>(constraint_objs[1]));
    objects9[1].Obj = cAEFace102;
    objects9[1].SubType = NXOpen::CAE::CaeSetObjectSubTypeNone;
    objects9[1].SubId = 0;
    setManagerFixed->SetTargetSetMembers(0, NXOpen::CAE::CaeSetGroupFilterTypeGeomFace, objects9);

    NXOpen::Unit* degree(dynamic_cast<NXOpen::Unit*>(workSimPart->UnitCollection()->FindObject("Degrees")));

    std::vector<NXOpen::Fields::FieldVariable*> indepVarArray2(0);
    fieldExpression1->EditFieldExpression("0", millimeter, indepVarArray2, false);

    propertyTable4->SetScalarFieldPropertyValue("DOF1", fieldExpression1);

    std::vector<NXOpen::Fields::FieldVariable*> indepVarArray3(0);
    fieldExpression2->EditFieldExpression("0", millimeter, indepVarArray3, false);

    propertyTable4->SetScalarFieldPropertyValue("DOF2", fieldExpression2);

    std::vector<NXOpen::Fields::FieldVariable*> indepVarArray4(0);
    fieldExpression3->EditFieldExpression("0", millimeter, indepVarArray4, false);

    propertyTable4->SetScalarFieldPropertyValue("DOF3", fieldExpression3);

    std::vector<NXOpen::Fields::FieldVariable*> indepVarArray5(0);
    fieldExpression4->EditFieldExpression("0", degree, indepVarArray5, false);

    propertyTable4->SetScalarFieldPropertyValue("DOF4", fieldExpression4);

    std::vector<NXOpen::Fields::FieldVariable*> indepVarArray6(0);
    fieldExpression5->EditFieldExpression("0", degree, indepVarArray6, false);

    propertyTable4->SetScalarFieldPropertyValue("DOF5", fieldExpression5);

    std::vector<NXOpen::Fields::FieldVariable*> indepVarArray7(0);
    fieldExpression6->EditFieldExpression("0", degree, indepVarArray7, false);

    propertyTable4->SetScalarFieldPropertyValue("DOF6", fieldExpression6);

    std::vector<NXOpen::NXString> propertyValue3(0);
    propertyTable4->SetTextPropertyValue("description", propertyValue3);

    fixedBuilder->SetDestinationFolder(nullNXOpen_CAE_SimLbcFolder);

    NXOpen::CAE::SimBC* fixed;
    fixed = fixedBuilder->CommitAddBc();

    fixedBuilder->Destroy();

    //time step table
    NXOpen::CAE::ModelingObjectPropertyTable* modelingObjectPropertyTableTimeTable;
    modelingObjectPropertyTableTimeTable = caePart1->ModelingObjectPropertyTables()->CreateModelingObjectPropertyTable("Time Step", "NX NASTRAN - Structural", "NX NASTRAN", "Time Step1", 4);
    NXOpen::CAE::PropertyTable* timeTable;
    timeTable = modelingObjectPropertyTableTimeTable->PropertyTable();
    timeTable->SetBaseScalarWithDataPropertyValue("Time Increment", "1.0", second);
    timeTable->SetIntegerPropertyValue("Number of Time Steps", 1);

    //Strategy Parameters
    NXOpen::CAE::ModelingObjectPropertyTable* modelingObjectPropertyTableStrategy;
    modelingObjectPropertyTableStrategy = caePart1->ModelingObjectPropertyTables()->CreateModelingObjectPropertyTable("Strategy Parameters", "NX NASTRAN - Structural", "NX NASTRAN", "Strategy Parameters1", 5);
    NXOpen::CAE::PropertyTable* StrategyTable;
    StrategyTable = modelingObjectPropertyTableStrategy->PropertyTable();

    StrategyTable->SetIntegerPropertyValue("MSTAB", 1);
    StrategyTable->SetIntegerPropertyValue("LSEARCH", 1);

    //добавление таблиц к simulation
    NXOpen::CAE::SimSolution* simSolution1(dynamic_cast<NXOpen::CAE::SimSolution*>(simSimulation1->ActiveSolution()));
    NXOpen::CAE::PropertyTable* propertyTable8;
    propertyTable8 = simSolution1->PropertyTable();
    std::vector<NXOpen::CAE::NamedPropertyTable*> propertyValue4(1);
    propertyValue4[0] = modelingObjectPropertyTableTimeTable;
    propertyTable8->SetNamedPropertyTableArrayPropertyValue("Time Step Intervals", propertyValue4);
    propertyTable8->SetNamedPropertyTablePropertyValue("Strategy Parameters", modelingObjectPropertyTableStrategy);
    propertyTable8->SetBooleanPropertyValue("Large Strains", true);
    propertyTable8->SetBooleanPropertyValue("Large Displacements", true);

    NXOpen::CAE::PropertyTable* nastranPropertiesTable;
    nastranPropertiesTable = simSolution1->SolverOptionsPropertyTable();

    nastranPropertiesTable->SetIntegerPropertyValue("parallel", 6);
    //sim_flag = true;
//  }
  //else theUI->NXMessageBox()->Show("CAD error", NXOpen::NXMessageBox::DialogTypeInformation, "Сначала постройте CAD-модель");


//  // UF_terminate();

//  //if (sim_flag && cad_flag) {
     // NXOpen::CAE::SimSolution* simSolution1(dynamic_cast<NXOpen::CAE::SimSolution*>(simSimulation1->ActiveSolution()));
      NXOpen::CAE::SimSolveManager* theCAESimSolveManager = NXOpen::CAE::SimSolveManager::GetSimSolveManager(theSession);
      std::vector<NXOpen::CAE::SimSolution*> psolutions1(1);
      psolutions1[0] = simSolution1;
      int numsolutionssolved1;
      int numsolutionsfailed1;
      int numsolutionsskipped1;

      

      theCAESimSolveManager->SolveChainOfSolutions(psolutions1, NXOpen::CAE::SimSolution::SolveOptionSolve, NXOpen::CAE::SimSolution::SetupCheckOptionCompleteCheckAndOutputErrors, NXOpen::CAE::SimSolution::SolveModeBackground, &numsolutionssolved1, &numsolutionsfailed1, &numsolutionsskipped1);
      
 // }

  ////////////////////////////////////////////////get results
      bool resFlag = true;
      NXOpen::CAE::SolutionResult* solutionResult1;
      NXOpen::CAE::SimResultReference* simResultReference1;
      while (resFlag) {
          Sleep(30000);
          simResultReference1 = (dynamic_cast<NXOpen::CAE::SimResultReference*>(simSolution1->Find("Structural")));
          try {

              solutionResult1 = theSession->ResultManager()->CreateReferenceResult(simResultReference1);
          }
          catch (...) {
              //  theUI->NXMessageBox()->Show("presAtNode", NXOpen::NXMessageBox::DialogTypeInformation, "fail");
              continue;
          }


          NXOpen::CAE::ResultParameters* resultParameters1;
          resultParameters1 = theSession->ResultManager()->CreateResultParameters();

          NXOpen::CAE::Loadcase* loadcase1(dynamic_cast<NXOpen::CAE::Loadcase*>(solutionResult1->Find("Loadcase[1]")));
          NXOpen::CAE::Iteration* iteration1(dynamic_cast<NXOpen::CAE::Iteration*>(loadcase1->Find("Iteration[1]")));
          NXOpen::CAE::ResultType* resultType1(dynamic_cast<NXOpen::CAE::ResultType*>(iteration1->Find("ResultType[[Contact Pressure][Nodal]]")));
          std::vector<NXOpen::CAE::Result::Section> section1;
          section1 = resultType1->GetSectionDefined();

          resultParameters1->SetDBScaling(0);

          NXOpen::CAE::SignalProcessingDBSettings* signalProcessingDBSettings1;
          signalProcessingDBSettings1 = resultParameters1->GetDbSettings();

          resultParameters1->SetGenericResultType(resultType1);

          resultParameters1->SetBeamSection(NXOpen::CAE::Result::BeamSection(-1));

          resultParameters1->SetShellSection(NXOpen::CAE::Result::ShellSection(-1));

          resultParameters1->SetResultComponent(NXOpen::CAE::Result::Component(-1));

          resultParameters1->SetCoordinateSystem(NXOpen::CAE::Result::CoordinateSystemAbsoluteRectangular);

          resultParameters1->SetSelectedCoordinateSystem(NXOpen::CAE::Result::CoordinateSystemSourceNone, -1);

          resultParameters1->SetRotationAxisOfAbsoluteCyndricalCSYS(NXOpen::CAE::Post::AxisymetricAxisNone);

          resultParameters1->SetBeamResultsInLocalCoordinateSystem(true);

          resultParameters1->SetShellResultsInProjectedCoordinateSystem(false);

          resultParameters1->MakeElementResult(false);

          resultParameters1->SetElementValueCriterion(NXOpen::CAE::Result::ElementValueCriterionAverage);

          resultParameters1->SetSpectrumScaling(NXOpen::CAE::SignalProcessingPlotData::ScalingTypeUnknown);

          resultParameters1->SetAcousticWeighting(NXOpen::CAE::SignalProcessingPlotData::AcousticalWeightingNone);

          NXOpen::CAE::Result::Averaging average1;
          average1.DoAveraging = false;
          average1.AverageAcrossPropertyIds = true;
          average1.AverageAcrossMaterialIds = true;
          average1.AverageAcrossElementTypes = true;
          average1.AverageAcrossFeatangle = true;
          average1.AverageAcrossAnglevalue = 45.0;
          average1.IncludeInternalElementContributions = true;
          resultParameters1->SetAveragingCriteria(average1);

          resultParameters1->SetComputationType(NXOpen::CAE::Result::ComputationTypeNone);

          resultParameters1->SetComputeOnVisible(false);

          resultParameters1->SetComplexCriterion(NXOpen::CAE::Result::ComplexAmplitude);

          resultParameters1->SetPhaseAngle(0.0);

          resultParameters1->SetSectionPlyLayer(0, 0, 1);

          resultParameters1->SetPlyID(0);

          resultParameters1->SetPlyLocation(NXOpen::CAE::Result::PlyLocationMiddle);

          resultParameters1->SetScale(1.0);

          NXOpen::Unit* unit1(dynamic_cast<NXOpen::Unit*>(workSimPart->UnitCollection()->FindObject("MilliMeter")));
          resultParameters1->SetUnit(unit1);

          resultParameters1->SetAbsoluteValue(false);

          resultParameters1->SetTensorComponentAbsoluteValue(NXOpen::CAE::Result::TensorDerivedAbsoluteDerivedComponent);

          resultParameters1->SetCalculateBeamStrResults(false);

          resultParameters1->SetBeamFillets(true);

          resultParameters1->SetBeamFilletRadius(0.5);

          resultParameters1->DisplayMidnodeValue(true);

          resultParameters1->SetIsReferenceNode(false);

          resultParameters1->SetReferenceNode(NULL);

          NXOpen::CAE::CyclicSymmetricParameters* cyclicSymmetricParameters1;
          cyclicSymmetricParameters1 = resultParameters1->GetCyclicSymmetricParameters();

          cyclicSymmetricParameters1->SetResultOption(NXOpen::CAE::CyclicSymmetricParameters::GetResultOnOriginalModel);

          cyclicSymmetricParameters1->SetOriginalResultOption(NXOpen::CAE::CyclicSymmetricParameters::OriginalResultBySector);

          cyclicSymmetricParameters1->SetSectCriteria(NXOpen::CAE::CyclicSymmetricParameters::SectorCriteriaIndex);

          cyclicSymmetricParameters1->SetSectorValue(NXOpen::CAE::CyclicSymmetricParameters::ValueMaximum);

          cyclicSymmetricParameters1->SetEnvValue(NXOpen::CAE::CyclicSymmetricParameters::EnvelopeValueAverage);

          cyclicSymmetricParameters1->SetSectorIndex(1);

          std::vector<int> sectors1(0);
          cyclicSymmetricParameters1->SetSectorIndices(sectors1);

          NXOpen::CAE::AxiSymmetricParameters* axiSymmetricParameters1;
          axiSymmetricParameters1 = resultParameters1->GetAxiSymmetricParameters();

          axiSymmetricParameters1->SetResultOption(NXOpen::CAE::AxiSymmetricParameters::GetResultOnOriginalModel);

          axiSymmetricParameters1->SetRotationAxis(NXOpen::CAE::AxiSymmetricParameters::AxisOfRotationXAxis);

          axiSymmetricParameters1->SetAxiOptions(NXOpen::CAE::AxiSymmetricParameters::OptionsAtRevolveAngle);

          axiSymmetricParameters1->SetEnvelopeVal(NXOpen::CAE::AxiSymmetricParameters::EnvValAverage);

          axiSymmetricParameters1->SetRevolveAngle(0.0);

          axiSymmetricParameters1->SetStartRevolveAngle(0.0);

          axiSymmetricParameters1->SetEndRevolveAngle(360.0);

          axiSymmetricParameters1->SetNumberOfSections(40);

          resultParameters1->SetProjectOnNodeNormal(false);

          int postviewId1;
          postviewId1 = theSession->Post()->CreatePostviewForResult(0, solutionResult1, false, resultParameters1);

          theSession->ResultManager()->DeleteResultParameters(resultParameters1);
          NXOpen::CAE::ResultAccess* resa = theSession->ResultManager()->CreateResultAccess(postviewId1);

          std::vector<int> nodes;
       
          for (int i = 0; i < nodesObject1.size(); i++) {
              nodes.push_back(nodesObject1[i]->Label());
          }
        

          double angPres = 0;
          std::vector<double> presAtNode = resa->AskNodalResult(nodes);
          int nodesNum = presAtNode.size();
          for (int i = 0; i < nodesNum; i++) {
              angPres += presAtNode[i];
          }
          angPres = angPres / nodesNum;
         
          theUI->NXMessageBox()->Show("presAtNode", NXOpen::NXMessageBox::DialogTypeInformation, std::to_string(angPres));
          resFlag = false;
      }
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