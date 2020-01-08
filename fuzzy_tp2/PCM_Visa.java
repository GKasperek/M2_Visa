  import ij.*;
  import ij.process.*;
  import ij.gui.*;
  import java.awt.*;
  import ij.plugin.frame.*;
  import ij.process.ImageProcessor.*;
  import ij.plugin.*;
  import java.awt.*;
  import java.awt.event.*;
  import java.awt.image.*;
  import java.io.*;
  import java.util.*;
  import java.lang.Math.*;
  import java.lang.Object.*;
  import java.lang.String.*;
  import java.awt.TextArea;
  import javax.swing.*;
  import javax.swing.event.*;
  import java.awt.Window.*;


  public class PCM_Visa implements PlugIn
  {

    class Vec{
      int[] data = new int[3];	//*pointeur sur les composantes*/
    }

    ////////////////////////////////////////////////////
    Random r = new Random();
    public int rand(int min, int max){
      return min + (int)(r.nextDouble()*(max-min));
    }


    ////////////////////////////////////////////////////////////////////////////////////////////
    public void run(String arg){

      // LES PARAMETRES
      ImageProcessor ip;
      ImageProcessor ipseg;
      ImageProcessor ipJ;
      ImagePlus imp;
      ImagePlus impseg;
      ImagePlus impJ;
      IJ.showMessage("Algorithme PCM","If ready, Press OK");
      ImagePlus cw;

      imp = WindowManager.getCurrentImage();
      ip = imp.getProcessor();

      int width = ip.getWidth();
      int height = ip.getHeight();


      impseg=NewImage.createImage("Image segment�e par FCM",width,height,1,24,0);
      ipseg = impseg.getProcessor();
      impseg.show();


      int nbclasses,nbpixels,iter;
      double stab,seuil,valeur_seuil;
      int i,j,k,l,imax,jmax,kmax;

      String demande =JOptionPane.showInputDialog("Nombre de classes : ");
      nbclasses =Integer.parseInt(demande);
      nbpixels = width * height; // taille de l'image en pixels

      demande =JOptionPane.showInputDialog("Valeur de m : ");
      double m =Double.parseDouble(demande);

      demande =JOptionPane.showInputDialog("Nombre iteration max : ");
      int itermax =Integer.parseInt(demande);

      demande =JOptionPane.showInputDialog("Valeur du seuil de stabilite : ");
      valeur_seuil =Double.parseDouble(demande);

      demande =JOptionPane.showInputDialog("Randomisation amelioree ? ");
      int valeur=Integer.parseInt(demande);

      // Centroides
      double c[][] = new double[nbclasses][3];

      // Centroides precedent
      double cprev[][] = new double[nbclasses][3];
      int cidx[] = new int[nbclasses];
      //double m;
      // Distance
      double Dmat[][] = new double[nbclasses][nbpixels];
      // Distance precedente
      double Dprev[][] = new double[nbclasses][nbpixels];

      //Degré d'appartenance
      double Umat[][] = new double[nbclasses][nbpixels];

      // Degré d'appartenance précédent
      double Uprev[][] = new double[nbclasses][nbpixels];

      // n grec penality term
      double N_grec[] = new double[nbclasses];

      double red[] = new double[nbpixels];
      double green[] = new double[nbpixels];
      double blue[] = new double[nbpixels];

      double[][] Pix = new double[nbpixels][3];

      int[] colorarray = new int[3];
      int[] init=new int[3];
      double figJ[]=new double[itermax];
      for(i=0;i<itermax;i++)
      {
        figJ[i]=0;
      }

      // R�cup�ration des donn�es images
      l = 0;
      for(i = 0; i < width; i++)
      {
        for(j = 0; j < height; j++)
        {
          ip.getPixel(i,j,colorarray);
          red[l] = (double)colorarray[0];
          green[l] =(double) colorarray[1];
          blue[l] = (double)colorarray[2];
          Pix[l][0] = red[l];
          Pix[l][1] = green[l];
          Pix[l][2] = blue[l];
          l++;
        }
      }
      ////////////////////////////////
      // FCM
      ///////////////////////////////

      imax = nbpixels;  // nombre de pixels dans l'image
      jmax = 3;  // nombre de composantes couleur
      kmax=nbclasses;
      double data[][] = new double[nbclasses][3];
      int[] fixe=new int[3];
      int xmin = 0;
      int xmax = width;
      int ymin = 0;
      int ymax = height;
      int rx, ry;
      int x,y;
      int epsilonx,epsilony;


      // Initialisation des centroedes (aleatoirement )

      for(i=0;i<nbclasses;i++){
        if(valeur==1){
          epsilonx=rand((int)(width/(i+2)),(int)(width/2));
          epsilony=rand((int)(height/(4)),(int)(height/2));
        }else{
          epsilonx=0;
          epsilony=0;
        }
        rx = rand(xmin+epsilonx, xmax-epsilonx);
        ry = rand(ymin+epsilony, ymax-epsilony);
        ip.getPixel(rx,ry,init);
        c[i][0] = init[0];
        c[i][1] =init[1];
        c[i][2] = init[2];
      }

      // Calcul de distance entre data et centroides
      for(l = 0; l < nbpixels; l++)
      {
        for(k = 0; k < kmax; k++)
        {
          double r2 = Math.pow(red[l] - c[k][0], 2);
          double g2 = Math.pow(green[l] - c[k][1], 2);
          double b2 = Math.pow(blue[l] - c[k][2], 2);
          Dmat[k][l] = r2 + g2 + b2;
        }
      }

      // Initialisation des degr�s d'appartenance
      //A COMPLETER
      for(i = 0;i <nbclasses; i++){
        for(j = 0; j < nbpixels; j++){
          double apparDegre_ij = 0.0;
          apparDegre_ij = this.appartenance_ij(i,j,m,c,Dmat);
          Umat[i][j] = apparDegre_ij;
        }
      }

      for(i = 0; i < nbclasses; i++){
        N_grec[i] = this.nGrec_i(i, Umat, Dmat, m);
      }
      ////////////////////////////////////////////////////////////
      // FIN INITIALISATION FCM
      ///////////////////////////////////////////////////////////


      /////////////////////////////////////////////////////////////
      // BOUCLE PRINCIPALE
      ////////////////////////////////////////////////////////////
      iter = 0;
      stab = 2;
      seuil = valeur_seuil;


      /////////////////// A COMPLETER ///////////////////////////////
      while ((iter < itermax) && (stab > seuil))
      {
        cprev = c;

        c = new double[c.length][3];
        for(i = 0; i < nbclasses; i++){
          c[i] = this.v_i(i,Pix,Umat,m);
        }
        // Update  the matrix of centroids
      /*  for(k = 0; k < nbclasses; k++){
          double redC   = 0.0;
          double greenC = 0.0;
          double blueC  = 0.0;

          double den = 0.0;*/

      /*    for(i = 0; i < nbpixels; i++){
            redC   = redC + (Math.pow(Umat[k][i], m) * red[i]);
            greenC = greenC + (Math.pow(Umat[k][i], m) * green[i]);
            blueC  = blueC + (Math.pow(Umat[k][i], m) * blue[i]);
            den += Math.pow(Umat[k][i], m);
          }
          if(den > 0){
            c[k][0] = redC/ den;
            c[k][1] = greenC / den;
            c[k][2] = blueC / den;
          }
        }*/

        // Compute Dmat, the matrix of distances (euclidian) with the centro�ds
        for(j = 0; j < nbpixels; j++){
          double[] xj = Pix[j];
          for(i =0; i<nbclasses; i++){
            double[] vi = c[i];
            double dij = this.mydist(vi,xj);

            Dmat[i][j] = dij;
          }
        }

        // Calcule la ùatrice des degrés d'appartenance_ij
        for(j = 0; j < nbpixels;j++){
          for(i = 0; i <nbclasses; i++){
            double appart = this.appartenance_ij(i,j,m,c,Dmat);
            Umat[i][j] = appart;
          }
        }

        // Calcul de ni
        for(i = 0; i < nbclasses; i++){
          N_grec[i] = this.nGrec_i(i,Umat,Dmat,m);
        }


        figJ[iter] = this.perf(Umat, Dmat, m, N_grec);

        if(iter > 0 && Math.abs(figJ[iter - 1]) <= Math.abs(figJ[iter])){
          break;
        }

        stab = Math.abs(figJ[iter] - stab);

        iter++;

    /*  for(k = 0; k < kmax; k++){
          for( i = 0; i < nbpixels; i++){
            Dprev[k][i] = Dmat[k][i];
            Uprev[k][i] = Umat[k][i];
          }
        }*/
        ////////////////////////////////////////////////////////

        // Affichage de l'image segment�e
        double[] mat_array=new double[nbclasses];
        l = 0;
        for(i=0;i<width;i++)
        {
          for(j = 0; j<height; j++)
          {
            for(k = 0; k<nbclasses; k++)
            {
              mat_array[k]=Umat[k][l];
            }
            int indice= IndiceMaxOfArray(mat_array,nbclasses) ;
            int array[] = new int[3];
            array[0] = (int)c[indice][0];
            array[1] = (int)c[indice][1];
            array[2] = (int)c[indice][2];
            ipseg.putPixel(i, j, array);
            l++;
          }
        }
        impseg.updateAndDraw();
        //////////////////////////////////
      }  // Fin boucle

      double[] xplot= new double[itermax];
      double[] yplot=new double[itermax];
      for(int w = 0; w < itermax; w++)
      {
        xplot[w]=(double)w;	yplot[w]=(double) figJ[w];
      }
      Plot plot = new Plot("Performance Index (PCM)","iterations","J(P) value",xplot,yplot);
      plot.setLineWidth(2);
      plot.setColor(Color.blue);
      plot.show();
    } // Fin FCM
    int indice;
    double min,max;

    private double perf(double[][] U, double[][] dist, double m,double[] N){
      double perfcm = 0;
      double penal = 0;

      for(int i = 0; i <U.length; i++){
        double total = 0;
        for(int j = 0; j < U[0].length; j++){
          perfcm += Math.pow(U[i][j], m) * dist[i][j];
          total += Math.pow(1 - U[i][j], m);
        }
        penal += N[i] * total;
      }
      perfcm += penal;
      return perfcm;
    }


    private double nGrec_i(int i, double[][] U, double[][] D, double m){
      double n;

      double num = 0;
      double denom = 0;

      for(int j = 0; j < U[0].length; j++){
        num += Math.pow(U[i][j], m)*D[i][j];
        denom += Math.pow(U[i][j], m);
      }

      n = num/denom;
      return n;
    }

    private double mydist(double[] vi, double[] xj){
      double r2 = Math.pow(xj[0] - vi[0], 2);
      double g2 = Math.pow(xj[1] - vi[1], 2);
  		double b2 = Math.pow(xj[2] - vi[2], 2);
  		return r2 + g2 + b2;
    }

    private double[] v_i(int i, double[][] X, double[][] U, double m){
      double[] v = new double[] {0,0,0};

      double xR = 0;
      double xG = 0;
      double xB = 0;

      double Usum = 0;

      for(int j = 0; j < X.length; j++){
        xR += Math.pow(U[i][j], m) * X[j][0];
        xG += Math.pow(U[i][j], m) * X[j][1];
        xB += Math.pow(U[i][j], m) * X[j][2];

        Usum += Math.pow(U[i][j], m);
      }
      v[0] = xR / Usum;
      v[1] = xG / Usum;
      v[2] = xB / Usum;

      return v;
    }

    private double appartenance_ij(int i, int j, double m, double[][] centro, double[][] distance){
      double appar = 0;

      if(distance[i][j] == 0){
        return 1;
      }

      for(int k = 0; k < centro.length; k++){
        if(distance[k][j] == 0){
          return 0;
        }
        appar = appar + Math.pow(distance[i][j]/distance[k][j], (2/m-1));
      }
      appar = Math.pow(appar, -1);
      return appar;
    }

    //Returns the maximum of the array
    public int  IndiceMaxOfArray(double[] array,int val)
    {
      max=0;
      for (int i=0; i<val; i++)
      {
        if (array[i]>max){
          max=array[i];
          indice=i;
        }
      }
      return indice;
    }

}
