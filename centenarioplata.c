#include <stdio.h>
#include <stdint.h>
#include <math.h>

#define M_PI 3.14159265358979323846
#define N_ITERATIONS 24

const double angles[] = {
    0.78539816339745, 0.46364760900081, 0.24497866312686, 0.12435499454676,
    0.06241880999596, 0.03123983343027, 0.01562372862048, 0.00781234106010,
    0.00390623013197, 0.00195312251648, 0.00097656218956, 0.00048828121119,
    0.00024414062015, 0.00012207031189, 0.00006103515617, 0.00003051757812,
    0.00001525878906, 0.00000762939453, 0.00000381469727, 0.00000190734863,
    0.00000095367432, 0.00000047683716, 0.00000023841858, 0.00000011920929,
    0.00000005960464, 0.00000002980232, 0.00000001490116, 0.00000000745058
};

const double kvalues[] = {
    0.70710678118655, 0.63245553203368, 0.61357199107790, 0.60883391251775,
    0.60764825625617, 0.60735177014130, 0.60727764409353, 0.60725911229889,
    0.60725447933256, 0.60725332108988, 0.60725303152913, 0.60725295913894,
    0.60725294104140, 0.60725293651701, 0.60725293538591, 0.60725293510314,
    0.60725293503245, 0.60725293501477, 0.60725293501035, 0.60725293500925,
    0.60725293500897, 0.60725293500890, 0.60725293500889, 0.60725293500888
};


struct ClarkeTransform{
    float alpha;
    float beta;
};

const struct ClarkeTransform Clarke(float a, float b, float c){
    return (const struct  ClarkeTransform)
    {
        .alpha = (2.0f/3.0f)*(a - b/2.0f - c/2.0f),
        .beta  = (2.0f/3.0f)*((sqrt(3.0f)/2.0f)*b - (sqrt(3.0f)/2.0f)*c)
    };
    
}

void cordic(double angle, int n, double *cosval, double *sinval){
    int i, ix, sigma;
    double kn, x, y, atn, t, theta = 0.0, pow2 = 1.0;
    int newsgn = (int)floor(angle / (2.0 * M_PI)) % 2 == 1 ? 1 : -1;
    

    if (angle < -M_PI/2.0 || angle > M_PI/2.0){
        if(angle < 0){
            cordic(angle+M_PI, n, &x, &y);
            
        } else {
            cordic(angle-M_PI, n, &x, &y);
        }
        *cosval = x * newsgn;
        *sinval = y * newsgn;
        return;
    }
    ix = n - 1;
    if(ix > 23) ix = 23;
    kn = kvalues[ix];
    x = 1;
    y =0 ;
    
    for (i = 0; i < n; ++i) {
        atn = angles[i];
        sigma = (theta < angle) ? 1 : -1;
        theta += sigma * atn;
        t = x;
        x -= sigma * y * pow2;
        y += sigma * t * pow2;
        pow2 /= 2.0;
    }
    *cosval = x * kn;
    *sinval = y * kn;    
}

struct sincosCalc
{
    float sin;
    float cos;
};

const struct sincosCalc trigCalc(double angle, int n){
    double coscalc, sincalc;
    
    cordic(angle, n, &coscalc, &sincalc);

    // Return sine and cosine as a const struct
    return (const struct sincosCalc){ .cos = coscalc, .sin = sincalc };
}


struct parkTransform
{
    float sd;
    float sq;
};

const struct parkTransform Park(float alph, float bet, float theta){
    
    struct sincosCalc result = trigCalc(theta, 16);
    
    return (const struct parkTransform)
    {.sd = alph*result.cos + bet*result.sin, 
     .sq =  -alph*result.sin+bet*result.cos};
}



struct PIDController {
    double Kp;        // Proportional gain
    double Ki;        // Integral gain
    double Kd;        // Derivative gain
    double prevError; // Previous error value
    double integral;  // Integral of error
};

// Function to initialize the PID controller
void PID_Init(struct PIDController *pid, double Kp, double Ki, double Kd) {
    pid->Kp = Kp;
    pid->Ki = Ki;
    pid->Kd = Kd;
    pid->prevError = 0.0;
    pid->integral = 0.0;
}

// Function to compute the PID output
double PID_Compute(struct PIDController *pid, double setpoint, double measurement, double dt) {
    
    
    // Calculate error
    double error = setpoint - measurement;

    // Proportional term
    double Pout = pid->Kp * error;

    // Integral term
    pid->integral += error * dt;
    double Iout = pid->Ki * pid->integral;

    // Derivative term
    double derivative = (error - pid->prevError) / dt;
    double Dout = pid->Kd * derivative;

    // Save current error for next iteration
    pid->prevError = error;

    // Calculate total output
    return Pout + Iout + Dout;
}

struct inversePark
{
    float iAlpha;
    float iBeta;
};

const struct inversePark invPark(float q, float d, float theta, double setPointq, double referenceq, double setPointd, double referenced){
    struct sincosCalc result = trigCalc(theta, 16);
    struct PIDController pidq;
    struct PIDController pidd;

    PID_Init(&pidq, 10, 5, 0);
    double fdback_q = PID_Compute(&pidq, setPointq, referenceq, 0.01);

    PID_Init(&pidd, 10, 5, 0);
    double fdback_d = PID_Compute(&pidd, setPointd, referenced, 0.01);

    return (const struct inversePark)
    {.iAlpha = d*result.cos - q*result.sin, 
     .iBeta =  d*result.sin + q*result.cos};
}

struct inverseClark
{
    double fdbackA;
    double fdbackB;
    double fdbackC;
};

const struct inverseClark invClark(double ia, double ib){
    return (const struct inverseClark)
    {
        .fdbackA = ia,
        .fdbackB = -(1.0f/2.0f)*ia + (sqrt(3.0f)/2.0f)*ib, 
        .fdbackC = -(1.0f/2.0f)*ia - (sqrt(3.0f)/2.0f)*ib
    };
    
}

int main(){
    
    const double amp =      11.0;
    const double freq =     50;  //freq in hz
    const double sampling = 1000; //samples p sec
    const double duration = 1; //sec
    const int numsamples  = (int)(sampling*duration);

    FILE *csvFile = fopen("output.csv", "w");
    if (csvFile == NULL) {
        perror("Error opening file");
        return 1;
    }

    // Write the header row to the CSV file
    fprintf(csvFile, "A, B, C, theta, alpha, beta, sd, sq, invA, invB\n");

    for (int i = 0; i < numsamples; i++) {
        double time = i / sampling;                     // Time in seconds
        const double input =  amp * sin(2.0f * M_PI * freq * time); // Sinusoidal input
        const double input1 = amp * sin(2.0f * M_PI * freq * time + M_PI/2.0f); // Sinusoidal input
        const double input2 = amp * sin(2.0f * M_PI * freq * time + M_PI/4.0f); // Sinusoidal input
        const double thetaGen = 1.0 * sin(2.0f * M_PI * 2.61 * time + M_PI/4.0f); // Sinusoidal input
        double angle2calc = 2 * M_PI * freq * time;
        //struct sincosCalc trigResult = trigCalc(angle2calc, 24);
        double spq, refq, spd, refd;

        // Test functions with the current input
        struct ClarkeTransform resultclarke = Clarke(input, input1, input2);
        struct parkTransform resultpark = Park(resultclarke.alpha, resultclarke.beta, thetaGen);
        struct inversePark resultInvPark = invPark(resultpark.sd, resultpark.sq, thetaGen, spq, refq, spd, refd);
        struct inverseClark resultClarke = invClark(resultInvPark.iAlpha, resultInvPark.iBeta);
        //fprintf(csvFile, "%f,%f,%f,%f,%f,%f\n", input, input1, input2, thetaGen, resultclarke.alpha, resultclarke.beta);
        fprintf(csvFile, "%f,%f,%f,%f,%f, %f,%f,%f,%f,%f, %f,%f,%f\n", input, input1, input2, thetaGen, resultclarke.alpha, 
            resultclarke.beta, resultpark.sd, resultpark.sq, resultInvPark.iAlpha, resultInvPark.iBeta, 
            resultClarke.fdbackA, resultClarke.fdbackB, resultClarke.fdbackC);
        //fprintf(csvFile, "%f,%f,%f,%f,%f,%f,%f,%f\n", input, input1, input2, calcalpha, calcbeta, calcd, calcq, thetaGen);
    }//10

    fclose(csvFile); // Close the file

    return 0;
}






