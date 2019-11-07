
#include "nelsonc.h"

#include <cstdio>
#include <ctime>


using namespace std;
double inz[40000];

nelsonc::nelsonc(int runs)
{
    n=runs;
    x0=0.0;
    //for i=1 to n ...
    rstart();
    //end for

    /* Testausgabe */
    cout<<"/******* Start: Testausgabe *******/\n"<<"\n";
    cout<<"        x0: "<<pos[0]<<"\n";
    cout<<"        g : "<<gauss(0.0,0.1);
    cout<<"        dt: "<<dt<<"\n";
    cout<<"\n\n/******* Ende: Testausgabe *******/\n"<<"\n";

}
bool ttt=true;
/* Startposition eines e- berechnen --> pos[0] */
void nelsonc::rstart()
{

    srand(time(NULL));


    /* numbers */

    double y = (rand()%997)/1000.0;
if(ttt)
{
    /*  Untere Grenze des Integrals */
    start[1]=5.0/(4.0*M_PI*a)*1.0/2.0*sqrt(M_PI*a)*exp(-s*s/a)*(exp(s*s/a)*(-2.0)-2.0);

    /* Integralwerte in Array ordnen   */
   for (int i = 0; i <= 20000; i++ )
   {
    r=-100.0+(i*0.005);
    inz[i]=5.0/(4.0*M_PI*a)*1.0/2.0*sqrt(M_PI*a)*exp(-s*s/a)*(exp(s*s/a)*(erf((r-s)/sqrt(a))+erf((r+s)/sqrt(a)))+2*erf(r/sqrt(a)))-start[1];

   }
   for (int i = 20001; i < 40000; i++ )
   {
    r=0.000+((i-20001)*0.005);
    inz[i]=5.0/(4.0*M_PI*a)*1.0/2.0*sqrt(M_PI*a)*exp(-s*s/a)*(exp(s*s/a)*(erf((r-s)/sqrt(a))+erf((r+s)/sqrt(a)))+2*erf(r/sqrt(a)))-start[1];

   }
ttt=false;
}
   /* y mit Integralwerten im Array vergleichen und daraus x0 bestimmen --> Startposition */
  if(y<0.498678)
  {
      for (int i = 0; i <= 20000; i++ )
      {
       if(y<inz[i]){break;}
       x0=-100.0+((double)i*0.005);
      }
  }
  else
  {
      for (int i = 20001; i < 40000; i++ )
      {
       if(y<inz[i]){break;}
       x0=0.0+(((double)i-20001)*0.005);
      }

  }

  pos[1]=0.0;
  pos[0]=x0;
  ti=0.0;


}

/* gaußverteilte Zufallszahl ausgeben */
double nelsonc::gauss(double mean,double stddev)
{
    random_device rd;
    mt19937 e2(rd());
    normal_distribution<> dist(mean, stddev);

    return  dist(e2);
}



void nelsonc::schritt()
{
    std::clock_t start;
    double duration;

       start = std::clock();
       cout<<"\nZeitmessung gestartet\n";


    ofstream myfile("pfad.txt"); // für die teilchenpfade
    double ruckstop;
    int histo[2001];
        int histol[2001];
            int histor[2001];
    double ddt[2001];
    for (int i = 0; i < 2001; i++ ){histo[i]=0;histol[i]=0;histor[i]=0;ddt[i]=0.0;}
    int maxx=10000;
    for (int i = 0; i <= 1000000; i++ )
    {
     rstart();

    while(pos[1]<600)
    {
       ruckstop=pos[1];
       pos[0]=pos[0]+(pos[0]*((ti-1.0)/(1.0+ti*ti))-((sinh(2.0*s*pos[0]/((a+a*ti*ti)))+sin((2.0*ti*s*pos[0]/(a+a*ti*ti))-(B*pos[1]*s/2.0)))/(cosh(2.0*s*pos[0]/(a+a*ti*ti))+cos((2.0*ti*s*pos[0]/(a+a*ti*ti))-(B*pos[1]*s/2.0))))*(s*(ti-1.0)/(1.0+ti*ti)))*dt+gauss(0.0,1.0)*sqrt(dt)*sqrt(2.0);
       pos[1]=pos[1]+(((ti-1.0)/(1.0+(ti)*(ti)))*pos[1]+8.0*a*(1.0+ti)/(1.0+(ti)*(ti))+(B*s*a/2.0)*(sinh((2.0*s*pos[0])/(a+a*ti*ti))-sin(((2.0*ti*s*pos[0])/(a+a*ti*ti))-(B*pos[1]*s/2.0)))/(cosh((2.0*s*pos[0])/(a+a*ti*ti))+cos(((2.0*s*pos[0]*ti)/(a+a*ti*ti))-(B*pos[1]*s/2.0))))*dt+gauss(0.0,1.0)*sqrt(dt)*sqrt(2.0);
       ti=ti+dt;
       if(pos[1]<ruckstop){pos[1]=ruckstop;}
       myfile<<pos[0]<<" "<<pos[1]<<"\n";// für die teilchenpfade
   }

   for (int i = 0; i < 2001; i++ )
   {
       if((-100.0+i/10.0)>pos[0])
       {
           if(x0<0.0){histol[i]++;}else{histor[i]++;}
           histo[i]++;ddt[i]=ddt[i]+ti;break;
       }
   }

    myfile<<"\n";// für die teilchenpfade
    if(maxx>10000)//zwischenergebnisse sichern
    {
        duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;


        cout<<i<<" Teilchen wurden simuliert. Gesamtzeit: "<<duration<<"s\n";maxx=0;

        ofstream myh("h.txt");
        for (int i = 0; i < 2001; i++ )
        {
         myh<<(-100.0+i/10.0)<<" "<<histo[i]<<" "<<ddt[i]/histo[i]<<" "<<histol[i]<<" "<<histor[i]<<"\n";
        }
        myh.close();

    }

    x0=0.0;pos[0]=0.0;pos[1]=0.0;
    ti=0.0;
    maxx++;


    }

     myfile.close();// für die teilchenpfade

     ofstream myh("h.txt");//ergebnis
     for (int i = 0; i < 2001; i++ )
     {
       myh<<(-100.0+i/10.0)<<" "<<histo[i]<<" "<<ddt[i]/histo[i]<<" "<<histol[i]<<" "<<histor[i]<<"\n";
     }
     myh.close();


}

