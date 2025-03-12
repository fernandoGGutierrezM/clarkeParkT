#include <stdio.h>
#include <stdint.h>
#include <math.h>

//sin and cos using taylor to aprox


// Factorials precomputed as constants
#define FACTORIAL_3 6      // 3!
#define FACTORIAL_5 120    // 5!
#define M_PI 3.14159265358979323846

float betainv, alphainv, zeta;
double alpha, beta, sd, sq, calcalpha, calcbeta, calcd, calcq;

// Function to calculate sine using Taylor series
float sin_taylor(float x) {
    // Normalize angle to [-π, π]
    while (x > 3.14159265359f) x -= 6.28318530718f;
    while (x < -3.14159265359f) x += 6.28318530718f;

    // Normalize angle to [-π/2, π/2] using symmetry
    int sign = 1;
    if (x > 1.57079632679f) {
        x = 3.14159265359f - x;
        sign = -1;
    } else if (x < -1.57079632679f) {
        x = -3.14159265359f - x;
        sign = -1;
    }

    // Compute Taylor series terms
    float x2 = x * x;                  // x^2
    float term1 = x;                   // First term: x
    float term3 = -(x * x2) / FACTORIAL_3;   // Second term: -(x^3)/3!
    float term5 = (x * x2 * x2) / FACTORIAL_5; // Third term: (x^5)/5!

    // Sum terms
    return sign * (term1 + term3 + term5);
}


// Function to calculate factorial
double factorial(int n) {
    double fact = 1.0;
    for (int i = 1; i <= n; i++) {
        fact *= i;
    }
    return fact;
}

// Function to calculate cosine using Taylor series
double cosine_taylor(double x, int terms) {
    double result = 0.0;

    // Loop through the series terms
    for (int n = 0; n < terms; n++) {
        // Calculate each term: (-1)^n * x^(2n) / (2n)!
        double term = ((n % 2 == 0 ? 1 : -1) * pow(x, 2 * n) / factorial(2 * n));
        result += term;
    }

    return result;
}


void clarkeT(double a, double b, double c){
    alpha = (2.0/3.0)*(a - b/2.0 - c/2.0);
    beta = (2.0/3.0)*((sqrt(3.0)/2.0)*b - (sqrt(3.0)/2.0)*c);
}

/*
void clarkeT(double a, double b, double c){
    alpha = 1;
    //alpha = (2/3)*(a - b/2 - c/2); 
    beta = (1/(3))*a + (2/(3))*b;
	//beta  = (2.0/3.0)*((sqrt(3.0)/2.0)*b-(sqrt(3.0)/2.0)*c);
}\
*/
/*

struct ReturnType Clarke(float a, float b, float c)
{
    const float alpha = (2/3)*(a - b/2 - c/2); 
    const float beta  = (2/3)*((sqrt(3.0f)/2)*b-(sqrt(3)/2)*c);

};
*/

void park(double alph, double be, double theta){
	sd = alph*cosine_taylor(theta, 5) + be*sin_taylor(theta);
	sq = -alph*sin_taylor(theta) + be*cosine_taylor(theta, 5);
}

void inversePark(float d, float q, float theta){
    alphainv = d*cosine_taylor(theta, 5) + q*sin_taylor(theta);
	betainv = -d*sin_taylor(theta) + q*cosine_taylor(theta, 5);
}



float PID(float ref, float ctrl, float kp, float ki, float kd, float dt){
    float error = ref - ctrl;
    //float response = 0;

    double a = (kp + kd/dt) * error;
    double b = (-kp + (ki*dt) - (2*kd/0.02)) * error;
    double c = (kd/0.02) * error;

    float response = a+b+c;

    return response;
}

int main(){

    //sin 2.3 = 0.74570521217
    //cos 2.3 = 0.66627602128
    //float calcsin = sin_taylor(2.3);
    //float calccos = cosine_taylor(2.3, 5);

    //printf("value for sin is: %f and cos: %f/", calcsin, calccos);

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
    fprintf(csvFile, "Input,Function1_Output,Function2_Output\n");

    for (int i = 0; i < numsamples; i++) {
        double time = i / sampling;                     // Time in seconds
        double input =  amp * sin(2 * M_PI * freq * time); // Sinusoidal input
        double input1 = amp * sin(2 * M_PI * freq * time + M_PI/2); // Sinusoidal input
        double input2 = amp * sin(2 * M_PI * freq * time + M_PI/4); // Sinusoidal input
        double thetaGen = 1.0 * sin(2 * M_PI * 2.61 * time + M_PI/4); // Sinusoidal input
        // Test functions with the current input
        
        clarkeT(input, input1, input2);
        park(alpha, beta, thetaGen);
        //double output1 = alpha;
        //float output2 = beta;
        //float response = PID(1.2, calcalpha, 10, 5, 0, 0.01);
        //inversePark(sd, sq, thetaGen);

        // Write input and outputs to the CSV file
        //fprintf(csvFile, "%f,%f,%f\n", input, output1, output2);
        //printf("alpha val: %f", calcalpha);
        fprintf(csvFile, "%f,%f,%f,%f,%f,%f,%f,%f\n", input, input1, input2, calcalpha, calcbeta, calcd, calcq, thetaGen);
    }//10

    fclose(csvFile); // Close the file

    /*
    clarkeT(2.1, 3.3, 1.4);
    park(alpha, beta, 1.9);

    float ra, rb;
    ra = PID(1.2, 0.7, 12 ,22, 3, 0.2);
    rb = PID(1.7, 1.4, 14, 24, 2, 0.2);

    printf("calculated sd val: %f sq: %f \n", ra, rb);
    */


    return 0;
}
