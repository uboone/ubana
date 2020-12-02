/*This function calculates a new event weight based on the flash chi2 x
It's intent is to fix data / MC discrpencies observed in the NuCC filter, which uses the flash chi2 to select events
For more details, see https://microboone-docdb.fnal.gov/cgi-bin/private/RetrieveFile?docid=32909&filename=082720.pdf&version=1*/
double reweightFlashChi2(double x, bool useNewFit=true){
   double p[9];
   //Old factors. Includes dirt and ext in fit
   const double oldFactors = {1.99384917e-01, -1.05022405e-10,  6.34916308e-05, -3.44705817e-03, 7.32059307e-02, -5.91696006e-01,  2.14667463e+00, -1.02545380e+00, 3.65915734e-01};
   //New factors. Excludes dirt and ext in fit
   const double newFactors = { 1.99926228e-01,  7.88770729e-10,  7.07517919e-05, -3.78293911e-03, 8.04431506e-02, -6.42532837e-01,  2.35485751e+00, -1.01116606e+00, 2.14473378e-01};
   if(useNewFit)
      p = newFactors;
   else:
      p = oldFactors;   
   
   double returnVal = exp(-p[0]*x)*(p[1]*pow(x, 6) + p[2]*pow(x,5) + p[3]*pow(x,4) + p[4]*pow(x,3) + p[5]*pow(x,2) + p[6]*x + p[7]) + p[8] ;
   
  
   //This function can go below 0.0 for very small values of chi2.
   //If that is the case, set the event weight to 1
   return (returnVal >0 ? returnVal : 1.0);
}   
