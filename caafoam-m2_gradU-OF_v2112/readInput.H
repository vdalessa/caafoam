const bool pressArtDiff =
    runTime.controlDict().lookupOrDefault("pressArtDiff", false);

const bool convArtDiff =
    runTime.controlDict().lookupOrDefault("convArtDiff", false);


scalar sigma =
    runTime.controlDict().lookupOrDefault("sigma",1.);

scalar kp =
    runTime.controlDict().lookupOrDefault("kp",0.);

scalar ku =
    runTime.controlDict().lookupOrDefault("ku",0.);

scalar ducLevelPress =
    runTime.controlDict().lookupOrDefault("ducLevelPress",1.5);

scalar ducLevelConv =
    runTime.controlDict().lookupOrDefault("ducLevelConv",1.5);

scalar epsilon =
    runTime.controlDict().lookupOrDefault("epsilon",0.);


   
