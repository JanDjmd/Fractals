/*
  JWildfire - an image and animation processor written in Java 
  Copyright (C) 1995-2011 Andreas Maschke
  This is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser 
  General Public License as published by the Free Software Foundation; either version 2.1 of the 
  License, or (at your option) any later version.
 
  This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
  even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  Lesser General Public License for more details.
  You should have received a copy of the GNU Lesser General Public License along with this software; 
  if not, write to the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
  02110-1301 USA, or see the FSF site: http://www.fsf.org.
*/
package org.jwildfire.create.tina.variation;
import odk.lang.FastMath;

import org.jwildfire.create.tina.base.XForm;
import org.jwildfire.create.tina.base.XYZPoint;
import org.jwildfire.base.Tools;
//import org.jwildfire.create.tina.base.Layer;
import static org.jwildfire.base.mathlib.Complex;
import static org.jwildfire.base.mathlib.MathLib.*;

    
  public class Triplex {   
  private double per_fix;
  public double re;
  public double im;
  public double Z;
  public double save_re;
  public double save_im;
  public double save_Z;
  

  // Constructors
  public Triplex() {
    re = 0.;
    im = 0.;
    Z  = 0.;
    save_re = 0.;
    save_im = 0.;
    save_Z  = 0.;
    per_fix = 0.;
  }

  public Triplex(double Rp) {
    re = Rp;
    im = 0.;
    Z  = 0.;
    save_re = 0.;
    save_im = 0.;
    save_Z  = 0.;
    per_fix = 0.;
  }

  public Triplex(int Rp) {
    re = (double) Rp;
    im = 0.;
    Z  = 0.;
    save_re = 0.;
    save_im = 0.;
    save_Z  = 0.;
    per_fix = 0.;
  }

  public Triplex(double Rp, double Ip, double Zp) {
    re = Rp;
    im = Ip;
    Z  = Zp;
    save_re = 0.;
    save_im = 0.;
    save_Z  = 0.;
    per_fix = 0.;
  }

  public Triplex(Triplex zz) {
    re = zz.re;
    im = zz.im;
    Z  = zz.Z;
    save_re = 0.;
    save_im = 0.;
    save_Z  = 0.;
    per_fix = 0.;
  }

  // Basic, utils etc
  public void One() {
    re = 1.0;
    im = 0.0;
    Z  = 0.;
    
  }

  public void ImOne() {
    re = 0.0;
    im = 1.0;
    Z  = 0.;
  }

  public void Zero() {
    re = 0.0;
    im = 0.0;
	Z  = 0.;
  }

  public void Copy(Triplex zz) {
    re = zz.re;
    im = zz.im;
  	Z  = 0.;
  }

  public void Flip() { // swap re <> im
    double r2 = im;
    double i2 = re;
    double Z2 = Z;
    re = r2;
    im = i2;
    Z  = Z2;
  }

  public void Conj() {
    Z = -Z;
  }

  public void Neg() {
    re = -re;
    im = -im;
    Z  = -Z;
  }

  public double Sig() {
    if (re == 0)
      return 0.;
    if (re > 0)
      return 1.;
    return -1.;
  }

  public double Sig2() { // avoids returning 0
    if (re >= 0)
      return 1.;
    return -1.;
  }

  public double Mag2() {
    return re * re + im * im + Z * Z;
  }

  public double Mag2eps() { // <-------------- CHANGE
    return re * re + im * im + Z * Z + 1e-20;
  }

  public double MagInv() {
    double M2 = this.Mag2();
    return (M2 < 1e-100 ? 1.0 : 1.0 / M2);
  }

  // Save / Restore utils
  public void Save() { // saves the current state
    save_re = re;
    save_im = im;
  	save_Z  = im;
  }

  public void Restore() { // revert to prev state, 
    //keeping into save re & im the current
    re = save_re;
    im = save_im;
  	Z  = save_Z;
  }

  public void Switch() { // revert to prev state, 
    // keeping into save re & im the current
    double r2 = save_re;
    double i2 = save_im;
    double Z2 = save_Z;
    save_re = re;
    save_im = im;
    save_Z  = Z;
    re = r2;
    im = i2;
  	Z  = Z2;
  }

  public void Keep(Triplex zz) { // saves zz
    save_re = zz.re;
    save_im = zz.im;
  	save_Z  = zz.Z;
  }

  public Triplex Recall() { // gives you what was saved
    return new Triplex(save_re, save_im, save_Z);
  }
  

  public void NextPow() {
    this.Mul(this.Recall());
  }

  // Arith

  public void Sqr() {
    double r2 = re * re - im * im;
    double i2 = 2 * re * im;
    double Z2 = 2 * re * Z;
    re = r2;
    im = i2;
  	Z = Z2;
  }

  public void Recip() {
    double mi = this.MagInv();
    re = re * mi;
    im = -im * mi;
  }

  public void Scale(double mul) {
    re = re * mul;
    im = im * mul;
  }

  public void Mul(Triplex zz) {
    if (zz.im == 0.0) {
      this.Scale(zz.re);
      return;
    }
    double r2 = re * zz.re - im * zz.im;
    double i2 = re * zz.im + im * zz.re;
    re = r2;
    im = i2;
  }

  
  // No Trig version
  public void TriMul(Triplex lc, Triplex z1) {
	double r1 = Math.sqrt((z1.re * z1.re) + (z1.im * z1.im));
	double r2 = Math.sqrt((lc.re * lc.re) + (lc.im * lc.im));
	double temp = z1.re;
	double a = 1.0 - ((z1.Z * lc.Z) / (r1 * r2));
	z1.re = a * (z1.re * lc.re - z1.im * lc.im);
	z1.im = a * (lc.re * z1.im + temp * lc.im);
	z1.Z = r2 * z1.Z + r1 * lc.Z;
	return; 
  	}
  
  // Trig version
	public void TrigMul(Triplex a, Triplex b) {
	double ra = Math.sqrt(a.re *a.re + a.im * a.im + a.Z * a.Z);
	//double phia = this.Atan(); //atan(a.im, a.re); // azimuth
	double phia = Math.atan2(a.im, a.re); // azimuth
	
	double thetaa = 0.0;
	 
	//thetaa = this.Asin(); // asin(a.Z/ra);
	thetaa = Math.sin(a.Z/ra);
	 
	double rb = Math.sqrt(b.re*b.re + b.im*b.im + b.Z*b.Z);
	//double phib = this.Atan(); //atan(b.im,b.re); // azimuth
	double phib = Math.atan2(b.im, b.re); // azimuth
	double thetab = 0.0;
	 
	//thetab = this.Asin(); // asin(b.Z/rb);
	thetab = Math.sin(b.Z/rb);
	 
	double r = ra*rb;
	double phi = phia + phib;
	double theta = thetaa + thetab;
	b.re = r*Math.cos(theta)*Math.cos(phi);
	b.im = r*Math.cos(theta)*Math.sin(phi);
	b.Z  = r*Math.sin(theta);
	 
	return; // new Triplex(r*cos(theta)*cos(phi), r*cos(theta)*sin(phi), r*sin(theta));
	}

  
  public void Div(Triplex zz) {
    // divides by zz
    // except if zz is 0
    double r2 = im * zz.im + re * zz.re;
    double i2 = im * zz.re - re * zz.im;
    double M2 = zz.MagInv();
    re = r2 * M2;
    im = i2 * M2;
  }

  public void DivR(Triplex zz) { // reverse
    double r2 = zz.im * im + zz.re * re;
    double i2 = zz.im * re - zz.re * im;
    double M2 = this.MagInv();
    re = r2 * M2;
    im = i2 * M2;
  }

  public void Add(Triplex zz) {
    re += zz.re;
    im += zz.im;
  }

  public void AMean(Triplex zz) {
    this.Add(zz);
    this.Scale(0.5);
  }


  public void RootMeanS(Triplex zz) {
    Triplex PF = new Triplex(zz);
    PF.Sqr();
    this.Sqr();
    this.Add(PF);
    this.Scale(0.5);
    this.Pow(0.5);
  }


  public void GMean(Triplex zz) {
    this.Mul(zz);
    this.Pow(0.5);
  }


  public void Heronian(Triplex zz) {
    // Heronian mean
    Triplex  HM = new Triplex(this);
    HM.GMean(zz);
    this.Add(zz);
    this.Add(HM);
    this.Scale(0.333333333333333333333333);

  }


  public void HMean(Triplex zz) {
    // not infinite if one term is zero. expansion by Maxima
    double p2 = (zz.re + re);
    double q2 = (zz.im + im);
    double D = 0.5 * (p2 * p2 + q2 * q2);
    if (D == 0) {
      this.Zero();
      return;
    }
    D = 1.0 / D;
    p2 = this.Mag2();
    q2 = zz.Mag2();
    if (p2 * q2 == 0) {
      this.Zero();
      return;
    }
    re = D * (re * q2 + zz.re * p2);
    im = D * (im * q2 + zz.im * p2);
    Z  = D * (Z  * q2 + zz.Z  * p2);
    
  }
  

  public void Sub(Triplex zz) {
    re -= zz.re;
    im -= zz.im;
    Z  -= zz.Z;
  }

  public void SubR(Triplex zz) {
    re = zz.re - re;
    im = zz.im - im;
    Z  = zz.Z  - Z;
    
  }

  public void Inc() {
    re += 1.0;
  }

  public void Dec() {
    re -= 1.0;
  }

  public void PerFix(double v) {
    // fix atan2() period, set v to a random integer
    per_fix = Math.PI * v;
  }

  public void Pow(double exp) {
    // some obvious cases come at first
    // instant evaluation for (-2, -1, -0.5, 0, 0.5, 1, 2)
    if (exp == 0.0) {
      this.One();
      return;
    }
	
    //double ex = MathLib.fabs(exp);
	double ex = fabs(exp);
    if (exp < 0) {
      this.Recip();
    }
	
	
    if (ex == 0.5) {
      this.Sqrt();
      return;
    }
	
	
    if (ex == 1.0) {
      return;
    }
	
	
    if (ex == 2.0) {
      this.Sqr();
      return;
    }
	
	
    // In general we need sin, cos etc
    Triplex PF = this.ToP();
    PF.re = Math.pow(PF.re, ex);
    PF.im = PF.im * ex;
    PF.Z = Z;
    this.Copy(PF.UnP());
	
  }


  // Trascendent functions (slower than others)

  public double Radius() {
    return Math.hypot(re, im);
  }

  public double Arg() {
    return (per_fix + Math.atan2(im, re));
  }

    public double ArgZ() {
    return (per_fix + Math.atan2(Z, re));
  }

  
  public Triplex ToP() {
    return new Triplex(this.Radius(), this.ArgZ(), Z);
  }

  public Triplex UnP() {
    return new Triplex(re * Math.cos(im), re * Math.sin(im), Z);
  }

  public void Norm() {
    this.Scale(Math.sqrt(this.MagInv()));
  }

  public void Exp() {
    re = Math.exp(re);
    this.Copy(this.UnP());
  }

  public void SinH() {
    double rr = 0.0;
    double ri = 0.0;
    double er = 1.0;
    re = Math.exp(re);
    er /= re;
    rr = 0.5 * (re - er);
    ri = rr + er;
    re = Math.cos(im) * rr;
    im = Math.sin(im) * ri;
  }

  public void Sin() {
    this.Flip();
    this.SinH();
    this.Flip(); // Should be this simple!
  }

  public void CosH() {
    double rr = 0.0;
    double ri = 0.0;
    double er = 1.0;
    re = Math.exp(re);
    er /= re;
    rr = 0.5 * (re - er);
    ri = rr + er;
    re = Math.cos(im) * ri;
    im = Math.sin(im) * rr;
  }

  public void Cos() {
    this.Flip();
    this.CosH();
    this.Flip(); // Should be this simple!
  }

  public void Sqrt() {
    // uses the exact formula and the expansion...
    // sqrt(a+i b) = sqrt(Radius + re) + i * sqrt(Radius - re)
    double Rad = this.Radius();
    double sb = (im < 0) ? -1 : 1;
    im = sb * Math.sqrt(0.5 * (Rad - re));
    re = Math.sqrt(0.5 * (Rad + re));
    // Re always gt 0. Primary value...
    if (per_fix < 0)
      this.Neg();
  }
 
  public void Log() {
    // I use Mag2eps to tame singularity
    Triplex L_eps = new Triplex(0.5 * Math.log(this.Mag2eps()), this.ArgZ(), this.ArgZ());
    this.Copy(L_eps);
  }


  public void LMean(Triplex zz) {
    // Logarithmic mean is given by (b-a)/log(b/a)
    Triplex dab = new Triplex(this);
    Triplex lab = new Triplex(this);
    dab.Sub(zz);
    lab.Div(zz);
    lab.Log();
    dab.Div(lab);
    this.Copy(dab);
  }


  public void AtanH() {
    Triplex D = new Triplex(this);
    D.Dec();
    D.Neg();
    this.Inc();
    this.Div(D);
    this.Log();
    this.Scale(0.5);
  }

  public void AsinH() { // slower than AtanH!
    Triplex D = new Triplex(this);
    D.Sqr();
    D.Inc();
    D.Pow(0.5);
    this.Add(D);
    this.Log();
  }

  public void AcosH() { // slower than Atanh!
    Triplex D = new Triplex(this);
    D.Sqr();
    D.Dec();
    D.Pow(0.5);
    this.Add(D);
    this.Log();
  }

  public void AcotH() {
    this.Recip();
    this.AtanH();
  }

  public void AsecH() {
    this.Recip();
    this.AsinH();
  }

  public void AcosecH() {
    this.Recip();
    this.AcosH();
  }

  public void Atan() { // the Flip() cheat works in my tests
    this.Flip();
    this.AtanH();
    this.Flip();
  }

  public void Asin() {
    // the Flip() cheat works in my tests
    this.Flip();
    this.AsinH();
    this.Flip();
  }

  public void Acos() {
    // not so suitable for fractals I think 
    // There is another one for this ... shift the Asin()
    this.Flip();
    this.AsinH();
    this.Flip(); // this is Asin()
    re = (Math.PI/2) - re;
    im = -im; // Acos = pi/2 - Asin
  }

  public void CPow(Triplex ex) {
    if (ex.im == 0.0) {
      this.Pow(ex.re);
      return;
    }
    this.Log();
    this.Mul(ex);
    this.Exp();
  }

}


public class SummerGardenFunc extends VariationFunc {
  private static final long serialVersionUID = 1L;
  
  private static final String PARAM_DIVIDE = "divide";
  private static final String PARAM_EXPAND = "expand"; 
  private static final String PARAM_MODE = "mode";    
  private static final String[] paramNames = { PARAM_DIVIDE, PARAM_EXPAND, PARAM_MODE};
  private double divide = 2.0;
  private double expand = 2.0;
  private int mode = 0;
  
  @Override
  public void transform(FlameTransformationContext pContext, XForm pXForm, XYZPoint pAffineTP, XYZPoint pVarTP, double pAmount) {
  // From Whittaker Courtney SummerGarden Mods Jan Dabo with Triplex Lib 
  // mode 0 default, mode 1 = swap sin and cos, mode 2 = square last z value.
  // default formula = sin(acos(z+2)/2) + cos(acos(z-2)/2)

       Triplex z = new Triplex(pAffineTP.x, pAffineTP.y, pAffineTP.z);
       Triplex zc = new Triplex(pAffineTP.x, pAffineTP.y, pAffineTP.z);       

       Triplex z5 = new Triplex(divide, 0, 0);
       Triplex z6 = new Triplex(expand, 0, 0);                
   
z.Add(z6);
zc.Sub(z6); 

z.Acos();
z.Div(z5);
//swap sin and cos depending on mode.
if (mode == 0 || mode == 2){
z.Sin();
}
else if (mode == 1){
z.Cos();
    }
    
zc.Acos();
zc.Div(z5);
//swap sin and cos depending on mode.
if(mode == 0 || mode == 2){
zc.Cos();
}
else if (mode == 1){
zc.Sin();    
    }

z.Add(zc);
      
//mirror horizontally or horizontally and vertically depending on mode
	if (mode == 2){
		z.Sqr();
	   if (pContext.random() < 0.5){                              
		   pVarTP.x += pAmount * z.re;  
		   pVarTP.y += pAmount * z.im;   
		   pVarTP.z += pAmount * z.Z;
		   
		}
		else {     
		   pVarTP.x += pAmount * -z.re; 
		   pVarTP.y += pAmount * -z.im;
			pVarTP.z += pAmount * z.Z;
		}
	}

	else{                                                                
	   if (pContext.random() < 0.25){                              
		   pVarTP.x += pAmount * z.re;  
		   pVarTP.y += pAmount * z.im;    
			pVarTP.z += pAmount * z.Z;
		}
		else if(pContext.random() > 0.25 && pContext.random() < 0.5){     
		   pVarTP.x += pAmount * -z.re; 
		   pVarTP.y += pAmount * -z.im;
			pVarTP.z += pAmount * -z.Z;
		} 
		else if(pContext.random() > 0.5 && pContext.random() < 0.75){     
		   pVarTP.x += pAmount * z.re; 
		   pVarTP.y += pAmount * -z.im;
			pVarTP.z += pAmount * -z.Z;
		} 
		else{     
		   pVarTP.x += pAmount * -z.re; 
		   pVarTP.y += pAmount * z.im;
			pVarTP.z += pAmount * z.Z;
		}                                       
	}




	if (pContext.isPreserveZCoordinate()) {
		pVarTP.z += pAmount * pAffineTP.z;
	}

  }
  
  @Override
  public String[] getParameterNames() {
    return paramNames;
  }

  @Override
  public Object[] getParameterValues() {
    return new Object[] {divide, expand, mode};
  }

  @Override
  public void setParameter(String pName, double pValue) {
    if (PARAM_DIVIDE.equalsIgnoreCase(pName))
      divide = pValue;
  
    else if (PARAM_EXPAND.equalsIgnoreCase(pName))
      expand = pValue;            
    else if (PARAM_MODE.equalsIgnoreCase(pName))
      mode = (int) Tools.limitValue(pValue, 0, 2);   
    else
      throw new IllegalArgumentException(pName);
  }

  @Override
  public String getName() {
    return "SummerGardenFunc";
  }

}