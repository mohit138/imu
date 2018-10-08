void alphabeta_init(Alphabeta *pyaw, float Alpha,float Beta,float THeta,float Rk ,float THeta_previous ,float THetaintegrate ,float THetaintegrate_previous){
 pyaw->alpha=Alpha;
 pyaw->beta=Beta;
 //pyaw->Theta=THeta;
 pyaw->rk=Rk;
 pyaw->Theta_previous=THeta_previous;
 //pyaw->Thetaintegrate=THetaintegrate;
 pyaw->Thetaintegrate_previous=THetaintegrate_previous;
}

void alphabeta_update(Alphabeta *pyaw,float measured_value){
    pyaw->Theta = pyaw->Theta_previous + ( pyaw->Thetaintegrate_previous * Time );
    pyaw->Thetaintegrate = pyaw->Thetaintegrate_previous;
    pyaw->rk = measured_value - pyaw->Theta;
    pyaw->Theta += pyaw->alpha * pyaw->rk;
    pyaw->Thetaintegrate += ( pyaw->beta * pyaw->rk )/Time;
    pyaw->Theta_previous = pyaw->Theta;
    pyaw->Thetaintegrate_previous = pyaw->Thetaintegrate;  
}
