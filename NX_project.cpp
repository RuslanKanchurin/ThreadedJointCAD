﻿#include "Joint.h"
using namespace NXOpen::BlockStyler;

NXOpen::Session* theSession;
NXOpen::Part* workPart;
NXOpen::Part* displayPart;



void ufusr(char* param, int* retcode, int paramLen){

        theSession = NXOpen::Session::GetSession();
        workPart = theSession->Parts()->Work();
        displayPart = theSession->Parts()->Display();

        Joint* joint1 = new Joint(theSession, workPart, displayPart);
        joint1->build();
        
        UF_terminate();
}
int ufusr_ask_unload(void)
{
    return (UF_UNLOAD_IMMEDIATELY);
}




