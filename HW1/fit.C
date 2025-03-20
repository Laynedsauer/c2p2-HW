
#include "fit.h"
#include <ap_int.h>
#include <ap_fixed.h> // Using fixed-point arithmetic for resource optimization
#include <hls_stream.h>
#include <ap_axi_sdata.h>
#include <iostream>  // For debug prints
#include <bitset>    // For binary output
#include "hls_math.h"
using namespace std;


void fit(int* tempxs, int* tempys, int* tempsigmas, int* templasts, int tempN) {
//#pragma HLS PIPELINE

	int chisq[6]={0};
	int sigmaa[6]={0};
	int sigmab[6]={0};
	int interb[6]={0};
	int b[6]={0};
	int a[6]={0};
	int event=0;
	int totalN[6]={0};
    int MAX = 500;
    int totalx;
    int totaltempN=0;
    int Sx[6]={0};
    int Sy[6]={0};
    int S[6]={0};
    int Stt[6]={0};
    int xevents[6][MAX];
    int yevents[6][MAX];
    int ti[6][MAX];

    int sigmayevents[6][MAX];
    int xvals[MAX];
    int yvals[MAX];
    int sigmayvals[MAX];
    int I = 0;
    int i = 0;
    int colIndex = 0; // Column index for current row


    for (i = 0; i <= tempN-1;) {

        int x = tempxs[i];
        int y = tempys[i];
        int y_unc = tempsigmas[i];
        int last = templasts[i];

        // Ensure that colIndex doesn't exceed MAX
        if (colIndex <= MAX) {
            // Assign the x, y, and y_unc values to the corresponding arrays
            xevents[I][colIndex] = x;
            yevents[I][colIndex] = y;
            sigmayevents[I][colIndex] = y_unc;
            Sx[I]=Sx[I]+x/(y_unc*y_unc);
            Sy[I]=Sy[I]+y/(y_unc*y_unc);
            S[I]=S[I]+1/(y_unc * y_unc);
            totalx=totalx+x;
            totaltempN=totaltempN+1;

        }

        // Debug output
        //std::cout << "xevents for I = " << I << " at colIndex = " << colIndex << " is " << xevents[I][colIndex] << std::endl;
        //std::cout << "yevents for I = " << I << " at colIndex = " << colIndex << " is " << yevents[I][colIndex] << std::endl;
        //std::cout << "sigmayevents for I = " << I << " at colIndex = " << colIndex << " is " << sigmayevents[I][colIndex] << std::endl;

        // If last == 1, increment I and reset colIndex for the next row
        if (last == 1) {
            std::cout << "I is " << I << std::endl;
            //std::cout << "first xevent from inside " << xevents[I][0] << std::endl;
            totalN[I]=colIndex;
            I++; // Move to the next row (I)
            colIndex = 0; // Reset column index for the new row
        }
        // Increment the column index if last != 1


   colIndex++;
   i++; }
//Doing Math


    for (int event = 0; event <= 5; event++) {
for (int point=0; point<= totalN[event]; point++){
	std::cout <<"    ti's " << event << '  '<< (1/sigmayevents[event][point])*(xevents[event][point]-(Sx[event]/S[event]))<< std::endl;

    	ti[event][point]=(1/sigmayevents[event][point]) * (xevents[event][point]-(Sx[event]/S[event]));
    	//Stt[event]=Stt[event]+ti[event][point];
    	std::cout <<"    ti's " << event << '  '<< ti[point][event]<< std::endl;

    }
}

    for (int event = 0; event <= 5; event++) {
    for (int point = 0; point<= totalN[event]; point++){
    	interb[event]=interb[event] + ti[event][point] * xevents[event][point]/sigmayevents[event][point];
        }
    b[event]=(1/Stt[event])*interb[event];
    a[event]=(Sy[event]-Sx[event]*b[event])/S[event];
    sigmaa[event]=(1/S[event])*(1+(Sx[event]*Sx[event])/(S[event]*Stt[event]));
    sigmab[event]=1/Stt[event];
    std::cout <<"b's " << b[event] << std::endl;
    std::cout <<"a's " << a[event] << std::endl;

    std::cout <<"sigmaa's " << sigmaa[event] << std::endl;

    std::cout <<"sigmab's " << sigmab[event] << std::endl;



    }

    for (int event = 0; event <= 5; event++) {
       for (int point = 0; point<=totalN[event]; point++){
       	chisq[event]=chisq[event]+(yevents[event][point]-a[event]-b[event]*xevents[event][point])/(sigmayevents[event][point]*(totalN[event]-2));
           }
       std::cout <<"chisq's " << chisq[event] << std::endl;
       }
	       std::cout <<"S[event]'s " << S[event] << std::endl;
	       std::cout <<"Stt[event]'s " << Stt[event] << std::endl;
    	   	   b[event]=(1/Stt[event])*interb[event];
    	       a[event]=(Sy[event]-Sx[event]*b[event])/S[event];
    	       sigmaa[event]=(1/S[event])*(1+(Sx[event]*Sx[event])/(S[event]*Stt[event]));
    	       sigmab[event]=1/Stt[event];
    	       std::cout <<"b's " << b[event] << std::endl;
    	       std::cout <<"a's " << a[event] << std::endl;

    	       std::cout <<"sigmaa's " << sigmaa[event] << std::endl;

    	       std::cout <<"sigmab's " << sigmab[event] << std::endl;



    }
