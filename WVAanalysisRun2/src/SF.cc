#include "math.h"

double SF_Photon (float eta, float pt) {

   double SF = 1.;

   if (pt >= 20 && pt < 30) {

      if (eta >= -2.5 && eta < -1.566) SF = 0.989799;
      else if (eta >= -1.4442 && eta < -1.) SF = 0.977669;
      else if (eta >= -1 && eta < 0) SF = 0.99272;
      else if (eta >= 0 && eta < 1.) SF = 0.984314;
      else if (eta >= 1. && eta < 1.4442) SF = 0.953749;
      else if (eta >= 1.566 && eta < 2.5) SF = 0.992262;

   } else if (pt >= 30 && pt < 40) {

      if (eta >= -2.5 && eta < -1.566) SF = 0.985586;
      else if (eta >= -1.4442 && eta < -1.) SF = 0.989285;
      else if (eta >= -1 && eta < 0) SF = 0.976581;
      else if (eta >= 0 && eta < 1.) SF = 0.980241;
      else if (eta >= 1. && eta < 1.4442) SF = 0.975098;
      else if (eta >= 1.566 && eta < 2.5) SF = 0.992262;

   } else if (pt >= 40 && pt < 50) {

      if (eta >= -2.5 && eta < -1.566) SF = 0.984684;
      else if (eta >= -1.4442 && eta < -1.) SF = 0.984485;
      else if (eta >= -1 && eta < 0) SF = 0.977651;
      else if (eta >= 0 && eta < 1.) SF = 0.978823;
      else if (eta >= 1. && eta < 1.4442) SF = 0.979511;
      else if (eta >= 1.566 && eta < 2.5) SF = 0.983761;

   } else if (pt >= 50 && pt < 200) {

      if (eta >= -2.5 && eta < -1.566) SF = 0.992006;
      else if (eta >= -1.4442 && eta < -1.) SF = 0.99123;
      else if (eta >= -1 && eta < 0) SF = 0.981936;
      else if (eta >= 0 && eta < 1.) SF = 0.975652;
      else if (eta >= 1. && eta < 1.4442) SF = 0.981635;
      else if (eta >= 1.566 && eta < 2.5) SF = 1.00227;

   }

   return (SF);

}

double SF_Elec (float eta, float pt) {

   double SF = 1.;

   if (pt >= 10 && pt < 20) {

      if (eta >= -2.5 && eta < -2.) SF = 1.221;
      else if (eta >= -2. && eta < -1.566) SF = 0.973;
      else if (eta >= -1.566 && eta < -1.4442) SF = 1.104;
      else if (eta >= -1.4442 && eta < -0.8) SF = 1.166;
      else if (eta >= -0.8 && eta < 0.) SF = 1.047;
      else if (eta >= 0. && eta < 0.8) SF = 1.035;
      else if (eta >= 0.8 && eta < 1.4442) SF = 1.129;
      else if (eta >= 1.4442 && eta < 1.566) SF = 1.311;
      else if (eta >= 1.566 && eta < 2.) SF = 1.078;
      else if (eta >= 2. && eta < 2.5) SF = 1.011;

   } else if (pt >= 20 && pt < 30) {

      if (eta >= -2.5 && eta < -2.) SF = 1.044;
      else if (eta >= -2. && eta < -1.566) SF = 0.996;
      else if (eta >= -1.566 && eta < -1.4442) SF = 1.193;
      else if (eta >= -1.4442 && eta < -0.8) SF = 0.972;
      else if (eta >= -0.8 && eta < 0.) SF = 0.894;
      else if (eta >= 0. && eta < 0.8) SF = 1.01;
      else if (eta >= 0.8 && eta < 1.4442) SF = 1.007;
      else if (eta >= 1.4442 && eta < 1.566) SF = 1.247;
      else if (eta >= 1.566 && eta < 2.) SF = 1.002;
      else if (eta >= 2. && eta < 2.5) SF = 1.032;

   } else if (pt >= 30 && pt < 40) {

      if (eta >= -2.5 && eta < -2.) SF = 0.986;
      else if (eta >= -2. && eta < -1.566) SF = 0.915;
      else if (eta >= -1.566 && eta < -1.4442) SF = 0.913;
      else if (eta >= -1.4442 && eta < -0.8) SF = 0.989;
      else if (eta >= -0.8 && eta < 0.) SF = 0.994;
      else if (eta >= 0. && eta < 0.8) SF = 0.965;
      else if (eta >= 0.8 && eta < 1.4442) SF = 0.971;
      else if (eta >= 1.4442 && eta < 1.566) SF = 0.945;
      else if (eta >= 1.566 && eta < 2.) SF = 0.935;
      else if (eta >= 2. && eta < 2.5) SF = 0.989;

   } else if (pt >= 40 && pt < 50) {

      if (eta >= -2.5 && eta < -2.) SF = 0.995;
      else if (eta >= -2. && eta < -1.566) SF = 1.001;
      else if (eta >= -1.566 && eta < -1.4442) SF = 1.02;
      else if (eta >= -1.4442 && eta < -0.8) SF = 0.983;
      else if (eta >= -0.8 && eta < 0.) SF = 0.971;
      else if (eta >= 0. && eta < 0.8) SF = 0.968;
      else if (eta >= 0.8 && eta < 1.4442) SF = 0.985;
      else if (eta >= 1.4442 && eta < 1.566) SF = 0.988;
      else if (eta >= 1.566 && eta < 2.) SF = 0.992;
      else if (eta >= 2. && eta < 2.5) SF = 0.991;

   } else if (pt >= 50 && pt < 200) {

      if (eta >= -2.5 && eta < -2.) SF = 0.998;
      else if (eta >= -2. && eta < -1.566) SF = 0.987;
      else if (eta >= -1.566 && eta < -1.4442) SF = 1.11;
      else if (eta >= -1.4442 && eta < -0.8) SF = 0.984;
      else if (eta >= -0.8 && eta < 0.) SF = 0.97;
      else if (eta >= 0. && eta < 0.8) SF = 0.989;
      else if (eta >= 0.8 && eta < 1.4442) SF = 0.98;
      else if (eta >= 1.4442 && eta < 1.566) SF = 1.096;
      else if (eta >= 1.566 && eta < 2.) SF = 0.985;
      else if (eta >= 2. && eta < 2.5) SF = 0.991;

   }

   return (SF);

}

double SF_Muon (float eta, float pt) {

   double SF = 1.;

   if (pt >= 20 && pt < 25) {

      if (fabs(eta) >= 0. && fabs(eta) < 0.9) SF = 0.975212;
      else if (fabs(eta) >= 0.9 && fabs(eta) < 1.2) SF = 0.97381;
      else if (fabs(eta) >= 1.2 && fabs(eta) < 2.1) SF = 0.998329;
      else if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) SF = 987784;

   } else if (pt >= 25 && pt < 30) {

      if (fabs(eta) >= 0. && fabs(eta) < 0.9) SF = 0.98483;
      else if (fabs(eta) >= 0.9 && fabs(eta) < 1.2) SF = 0.978645;
      else if (fabs(eta) >= 1.2 && fabs(eta) < 2.1) SF = 0.990546;
      else if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) SF = 0.980255;

   } else if (pt >= 30 && pt < 40) {

      if (fabs(eta) >= 0. && fabs(eta) < 0.9) SF = 0.986179;
      else if (fabs(eta) >= 0.9 && fabs(eta) < 1.2) SF = 0.979893;
      else if (fabs(eta) >= 1.2 && fabs(eta) < 2.1) SF = 0.992367;
      else if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) SF = 0.978505;

   } else if (pt >= 40 && pt < 50) {

      if (fabs(eta) >= 0. && fabs(eta) < 0.9) SF = 0.987443;
      else if (fabs(eta) >= 0.9 && fabs(eta) < 1.2) SF = 0.980234;
      else if (fabs(eta) >= 1.2 && fabs(eta) < 2.1) SF = 0.992763;
      else if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) SF = 0.977854;

   } else if (pt >= 50 && pt < 60) {

      if (fabs(eta) >= 0. && fabs(eta) < 0.9) SF = 0.983429;
      else if (fabs(eta) >= 0.9 && fabs(eta) < 1.2) SF = 0.97733;
      else if (fabs(eta) >= 1.2 && fabs(eta) < 2.1) SF = 0.988632;
      else if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) SF = 0.965441;

   } else if (pt >= 60 && pt < 120) {

      if (fabs(eta) >= 0. && fabs(eta) < 0.9) SF = 0.986318;
      else if (fabs(eta) >= 0.9 && fabs(eta) < 1.2) SF = 0.979523;
      else if (fabs(eta) >= 1.2 && fabs(eta) < 2.1) SF = 0.995045;
      else if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) SF = 0.968962;

   }

   return (SF);

}
