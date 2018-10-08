#include <WSWire.h>
#include <TimerOne.h>
#define IMU_V5
#define Time 0.1
// Uncomment the following line to use a MinIMU-9 v5 or AltIMU-10 v5. Leave commented for older IMUs (up through v4).

// Uncomment the below line to use this axis definition:
   // X axis pointing forward
   // Y axis pointing to the right
   // and Z axis pointing down.
// Positive pitch : nose up`
// Positive roll : right wing down
// Positive yaw : clockwise
int SENSOR_SIGN[9] = {1,1,1,-1,-1,-1,1,1,1}; //Correct directions x,y,z - gyro, accelerometer, magnetometer
// Uncomment the below line to use this axis definition:
   // X axis pointing forward
   // Y axis pointing to the left
   // and Z axis pointing up.
// Positive pitch : nose down
// Positive roll : right wing down
// Positive yaw : counterclockwise
//int SENSOR_SIGN[9] = {1,-1,-1,-1,1,1,1,-1,-1}; //Correct directions x,y,z - gyro, accelerometer, magnetometer

// tested with Arduino Uno with ATmega328 and Arduino Duemilanove with ATMega168
// accelerometer: 8 g sensitivity
// 3.9 mg/digit; 1 g = 256
#define GRAVITY 256  //this equivalent to 1G in the raw data coming from the accelerometer

#define ToRad(x) ((x)*0.01745329252)  // *pi/180
#define ToDeg(x) ((x)*57.2957795131)  // *180/pi

// gyro: 2000 dps full scale
// 70 mdps/digit; 1 dps = 0.07
#define Gyro_Gain_X 0.07 //X axis Gyro gain
#define Gyro_Gain_Y 0.07 //Y axis Gyro gain
#define Gyro_Gain_Z 0.07 //Z axis Gyro gain
#define Gyro_Scaled_X(x) ((x)*ToRad(Gyro_Gain_X)) //Return the scaled ADC raw data of the gyro in radians for second
#define Gyro_Scaled_Y(x) ((x)*ToRad(Gyro_Gain_Y)) //Return the scaled ADC raw data of the gyro in radians for second
#define Gyro_Scaled_Z(x) ((x)*ToRad(Gyro_Gain_Z)) //Return the scaled ADC raw data of the gyro in radians for second

// LSM303/LIS3MDL magnetometer calibration constants; use the Calibrate example from
// the Pololu LSM303 or LIS3MDL library to find the right values for your board
#define M_X_MIN -1198
#define M_Y_MIN -5962
#define M_Z_MIN +2286
#define M_X_MAX +5074
#define M_Y_MAX -438
#define M_Z_MAX +7739

#define Kp_ROLLPITCH 0.02
#define Ki_ROLLPITCH 0.00002
#define Kp_YAW 1.2
#define Ki_YAW 0.00002

/*For debugging purposes*/
//OUTPUTMODE=1 will print the corrected data,
//OUTPUTMODE=0 will print uncorrected data of the gyros (with drift)
#define OUTPUTMODE 1

#define PRINT_DCM 0     //Will print the whole direction cosine matrix
#define PRINT_ANALOGS 0 //Will print the analog raw data
#define PRINT_EULER 1   //Will print the Euler angles Roll, Pitch and Yaw

#define STATUS_LED 13

float G_Dt=0.02;    // Integration time (DCM algorithm)  We will run the integration loop at 50Hz if possible

long timer=0;   //general purpuse timer
long timer_old;
long timer24=0; //Second timer used to print values
int AN[6]; //array that stores the gyro and accelerometer data
int AN_OFFSET[6]={0,0,0,0,0,0}; //Array that stores the Offset of the sensors

int gyro_x;
int gyro_y;
int gyro_z;
int accel_x;
int accel_y;
int accel_z;
int magnetom_x;
int magnetom_y;
int magnetom_z;
float c_magnetom_x;
float c_magnetom_y;
float c_magnetom_z;
float MAG_Heading;

float Accel_Vector[3]= {0,0,0}; //Store the acceleration in a vector
float Gyro_Vector[3]= {0,0,0};//Store the gyros turn rate in a vector
float Omega_Vector[3]= {0,0,0}; //Corrected Gyro_Vector data
float Omega_P[3]= {0,0,0};//Omega Proportional correction
float Omega_I[3]= {0,0,0};//Omega Integrator
float Omega[3]= {0,0,0};

// Euler angles
float roll;
float pitch;
float yaw;

float errorRollPitch[3]= {0,0,0};
float errorYaw[3]= {0,0,0};

unsigned int counter=0;
byte gyro_sat=0;

float DCM_Matrix[3][3]= {  {1,0,0}  ,  {0,1,0}  , {0,0,1}  };
float Update_Matrix[3][3]={{0,1,2},{3,4,5},{6,7,8}}; //Gyros here


float Temporary_Matrix[3][3]={  {0,0,0}  ,  {0,0,0}  , {0,0,0}  };


#define KALMAN

#define MotorLeftDir1 A0
#define MotorLeftDir2 A2
#define MotorLeftPWM 8
#define MotorRightDir1 A3
#define MotorRightDir2 A7
#define MotorRightPWM 3
#define Kp 6
#define Kd 0.03
#define Ki 0
#define BasePWM_L 191
#define BasePWM_R 159
long int timer_begin;
float const_corr,angle,Difference_Error,Inte_Error,Prev_Error=0,measured_value;
int ref=0,PWM,LeftPWM,RightPWM;
float RightPWMIdeal;
float yaw_angle;
float timelag_check=0;


struct kalman
{
  float q;  //process noise covariance
  float r;  //measurement noise covariance
  float x;  //estimated value
  float p;  //error in estimate
  float k;  //kalman gain
}; 
typedef kalman Kalman;
Kalman raw;
Kalman *kf=&raw;

struct alphabeta{
  float alpha;
  float beta;
  float Theta;
  float rk; // prediction erroe
  float Theta_previous;
  float Thetaintegrate;
  float Thetaintegrate_previous;
};

typedef alphabeta Alphabeta;
Alphabeta yaw_ab;
Alphabeta *pyaw=&yaw_ab;


void update_kalman(Kalman *kf,float measured_value);

void kalman_init(Kalman *kf,float qq,float rr,float pp,float initial_estimate);

//void kalman_init(Kalman *praw,float qq,float rr,float pp,float initial_estimate)
//void update_kalman(Kalman *praw,measured_value)



int PID(float Error);


void setupIMU(void)
{
  Serial.begin(115200);
  pinMode (STATUS_LED,OUTPUT);  // Status LED

  I2C_Init();

  Serial.println("Pololu MinIMU-9 + Arduino AHRS");

  digitalWrite(STATUS_LED,LOW);
  delay(1500);

  Accel_Init();
  Compass_Init();
  Gyro_Init();

  delay(20);

  for(int i=0;i<32;i++)    // We take some readings...
    {
    Read_Gyro();
    Read_Accel();
    for(int y=0; y<6; y++)   // Cumulate values
      AN_OFFSET[y] += AN[y];
    delay(20);
    }

  for(int y=0; y<6; y++)  
    AN_OFFSET[y] = AN_OFFSET[y]/32;

  AN_OFFSET[5]-=GRAVITY*SENSOR_SIGN[5];

  //Serial.println("Offset:");
  for(int y=0; y<6; y++)
    Serial.println(AN_OFFSET[y]);

  delay(2000);
  digitalWrite(STATUS_LED,HIGH);

  timer=millis();
  delay(20);
  counter=0;
}

void getValues(void)
{
  
    counter++;
    timer_old = timer;
    timer=millis();
    if (timer>timer_old)
    {
      G_Dt = (timer-timer_old)/1000.0;    // Real time of loop run. We use this on the DCM algorithm (gyro integration time)
      if (G_Dt > 0.2)
        G_Dt = 0; // ignore integration times over 200 ms
    }
    else
      G_Dt = 0;



    // *** DCM algorithm
    // Data adquisition
    Read_Gyro();   // This read gyro data
    Read_Accel();     // Read I2C accelerometer

    if (counter > 5)  // Read compass data at 10Hz... (5 loop runs)
    {
      counter=0;
      Read_Compass();    // Read I2C magnetometer
      Compass_Heading(); // Calculate magnetic heading
    }

    // Calculations...
    Matrix_update();
    Normalize();
    Drift_correction();
    Euler_angles();
    // ***
    

   printdata();
}

int PID(float Error)
{
  float Final_Error;
  Difference_Error=Error-Prev_Error;
  Inte_Error+=Error;
  Final_Error=Kp*(Error)+Kd*(Difference_Error)+Ki*(Inte_Error);
  Prev_Error=Error;
  return (int)Final_Error;
}
char received[3];
int const_corr_rec;
void setup() {
  Serial.begin(9600);
  pinMode(MotorLeftDir1,OUTPUT);
  pinMode(MotorRightDir1,OUTPUT);
  pinMode(MotorLeftDir2,OUTPUT);
  pinMode(MotorRightDir2,OUTPUT);
  pinMode(MotorLeftPWM,OUTPUT);
  pinMode(MotorRightPWM,OUTPUT);
//  Serial.println("enter");

  
  kalman_init(kf,50,50,50,50);
  alphabeta_init(pyaw,0.99,0.005,1,1,0,0,0);

    
  setupIMU();

  //  Serial.println("exit");

  timer_begin=millis();
  while(1)
  {
    //Serial.println("atak gya");
    getValues();
    // printdata();
    // Serial.println(MAG_Heading*180/3.14);
    
    if(((ToDeg(roll)<1.00)&&(ToDeg(roll)>-1.00))&&((ToDeg(pitch)<1.00)&&(ToDeg(pitch)>-1.00))&&((millis()-timer_begin)>13000))
    {
       const_corr=ToDeg(yaw);
     // const_corr=MAG_Heading*180/3.14;
      
       break;
    }
    
    
  }
  
  //delay(1000);
 // const_corr=60;
  //Serial.println("working");
  digitalWrite(MotorLeftDir1,HIGH);
  digitalWrite(MotorLeftDir2,  LOW);
  digitalWrite(MotorRightDir1,HIGH);
  digitalWrite(MotorRightDir2,LOW);  
  
}

/*
int i=2;
bool stringComplete = false;
*/
void loop() {
/*
   receiveSerial();  */  
   if((millis()-timer)>=10)  // Main loop runs at 50Hz
  {
   getValues();
  // const_corr_rec = (received[2]-'0')*100 + (received[1]-'0')*10 + (received[0]-'0');
 //  const_corr=const_corr_rec-180;
  
 update_kalman(kf,ToDeg(yaw));   
  if(millis()-timelag_check>=3)
  {
      
   timelag_check=millis();
  alphabeta_update(pyaw,ToDeg(yaw));
  } 
 //         TEST FILTERS 
//    Serial.print(ToDeg(yaw));
//    L̥
//    Serial.print("     "+ String(kf->x));
//    
//    Serial.println("     "+ String(pyaw->Theta));
//    

  angle=(kf->x)-const_corr;

//  Serial.print(const_corr);
//  Serial.print(" "); 
//
//  Serial.println(kf->x);
//  
 
  if(angle<(-180))
  {
    angle=360-angle;
  }
  if(angle>(180))
  {
    angle=angle-360;
  } 
   PWM=PID(angle);
  // Serial.println(angle);
   LeftPWM=BasePWM_L-PWM;
   RightPWMIdeal=BasePWM_R+PWM;

    if(LeftPWM>255) LeftPWM=255;
    else if(LeftPWM<0) LeftPWM=0;

    if(RightPWMIdeal>255) RightPWMIdeal=255;
    else if(RightPWMIdeal<0) RightPWMIdeal=0;

  //  Serial.println(LeftPWM);
   // Serial.println(RightPWMIdeal);
    
   // RightPWMIdeal*=(0.8235);
    RightPWM=(int)RightPWMIdeal;
  // Serial.println(LeftPWM);
   // Serial.println(RightPWM);
   analogWrite(MotorLeftPWM,LeftPWM);
   analogWrite(MotorRightPWM,RightPWM);
  
  } 
 
}
/*
void receiveSerial()
{ while(Serial.available())
   {
    if(!stringComplete)
    { received[i]=Serial.read();
   //  Serial.println("Rec   "+String(received[i]));
     i--;
    }
     if(i<0){ i = 2;
     stringComplete = true;}
     
    }
    if(stringComplete)
    {
      const_corr_rec = (received[2]-'0')*100 + (received[1]-'0')*10 + (received[0]-'0');
      stringComplete = false;
    }
}*/
