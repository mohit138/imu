
float qq,rr,pp,initial_estimate;


void kalman_init(Kalman *kf,float qq,float rr,float pp,float initial_estimate){
  kf->q=qq;                          // kf kalman struct
                                     // x is estimate of sensor reading or input for the filter
  kf->r=rr;                           
  kf->p=pp;
  kf->x=initial_estimate;
}

void update_kalman(Kalman *kf,float measured_value){

  kf->p=kf->p+kf->q;                        //r is the error in measurement
  kf->k=kf->p/(kf->p + kf->r);              // k is the kalman gain
  kf->x=kf->x+kf->k*(measured_value-kf->x); //p is the error in estimate
  kf->p=(1-kf->k)*kf->p;                    //q process noise covariance
}

