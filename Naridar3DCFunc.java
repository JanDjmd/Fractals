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

import java.util.Random;
import org.jwildfire.create.tina.base.Layer;
import org.jwildfire.create.tina.base.XForm;
import org.jwildfire.create.tina.base.XYZPoint;
import org.jwildfire.create.tina.variation.FlameTransformationContext;
import org.jwildfire.base.Tools;
import static org.jwildfire.base.mathlib.MathLib.*;


public class Naridar3DCFunc extends VariationFunc {
// By Edgar Malinovsky Mods Jan Dabo for Jwildfire with Scale normally Scale=2
  private static final long serialVersionUID = 1L;
  
  private static final String PARAM_FoldX = "FoldX";
  private static final String PARAM_FoldY = "FoldY";
  private static final String PARAM_FoldZ = "FoldZ";
  private static final String PARAM_FoldW = "FoldW";
  private static final String PARAM_Scale = "Scale";

  private static final String PARAM_COLOR_MODE = "color mode 1-4"; 
  private static final String PARAM_COLOR_SCALE = "color scale";
	
  private static final String[] paramNames = {PARAM_FoldX, PARAM_FoldY, PARAM_FoldZ, PARAM_FoldW, PARAM_Scale, PARAM_COLOR_MODE, PARAM_COLOR_SCALE};
  
	private double FoldY = 1.0;
	private double FoldW = 1.0;
	private double FoldX = 1.0;
	private double FoldZ = 0.0;
	private double Scale = 2.0;
			
  	private int color_mode = 1;
  	private double color_scale = 1.0;


  private static double EPS = 1.0e-20;

	@Override
	public void transform(FlameTransformationContext pContext, XForm pXForm, XYZPoint pAffineTP, XYZPoint pVarTP, double pAmount) {
  
	// Parameters
		double xx=pAffineTP.x;
		double yy=pAffineTP.y;
		double zz=pAffineTP.z;

		double x = pAffineTP.x;
		double y = pAffineTP.y;
		double z = pAffineTP.z;
		double w = 1.0;
  
  // CODE GOES HERE

      	double x1, y1, z1, w1, x2, y2, z2, w2;
		double contrscale;
		double cx, cy, cz, cw;

		cx = pAffineTP.x;
		cy = pAffineTP.y;
		cz = pAffineTP.z;
		cw = 1.0;

		//=====================================
		///folds of different lenght. 
		x = cx +  Math.abs(cx-FoldX) - Math.abs(cx+FoldX);
		y = cy +  Math.abs(cy-FoldY) - Math.abs(cy+FoldY);
		z = cz +  Math.abs(cz-FoldZ) - Math.abs(cz+FoldZ);
		w = cw +  Math.abs(cw-FoldW) - Math.abs(cw+FoldW);	
		

//========================================
	x1 = x* Scale;
	y1 = y* Scale;
	z1 = z* Scale;
	w1 = w* Scale;


//===============================
//Now goes algebra. Calculates two formulas



// hyper tricorn - hypercomplex number tricorn.
	x2 = x1*x1 - y1*y1 - z1*z1 + w1*w1;
	y2 = Scale* z1*w1 - Scale*x1*y1;
	z2 = Scale*x1*z1 - Scale* y1*w1;
	w2 = Scale*x1*w1 + Scale* y1*z1;

//==================================================================

//quaternion pow2
	x = x1*x1 - y1*y1 - z1*z1 - w1*w1;
	y = Scale*x1*y1;
	z = Scale*x1*z1;
	w = Scale*w1*x1;


//=============================================================
//then choses larger value of the two

	if (x2>x) 
	x = x2;

	if (y2>y) 
	y = y2;

	if (z2>z) 
	z = z2;

	if (w2>w) 
	w = w2;

//========================================
//scaled C (z=z^2+c) generates strechless fractal
	contrscale = 1.00 /Scale;
	x = x + contrscale* pAffineTP.x; 
	y = y + contrscale* pAffineTP.y; 
	z = z + contrscale* pAffineTP.z; 
	

//
		pVarTP.x = pAmount * x;
		pVarTP.y = pAmount * y;
		pVarTP.z = pAmount * z;
		
		
		double mag1 = (sqrt(xx*xx+yy*yy+zz*zz)/2); //*color_scale;
       
		double mag2 = mag1 / color_scale;

		if (color_mode == 1){
		 pVarTP.color = mag2;
		}
		else if (color_mode == 2){
		 pVarTP.color += mag2;             
		}
		else if (color_mode == 3){
		 pVarTP.color = sin(sqrt(xx*xx+yy*yy+zz*zz)+mag2*2);             
		}       
		else if (color_mode == 4){
		 pVarTP.color += sin(sqrt(xx*xx+yy*yy+zz*zz)+mag2*2);             
		}  
		else if (color_mode == 5){
		 pVarTP.color = log(sqrt(xx*xx+yy*yy+zz*zz)+mag2*2);             
		}       
		else if (color_mode == 6){
		 pVarTP.color += log(sqrt(xx*xx+yy*yy+zz*zz)+mag2*2);             
		}  

  
  }
  
    @Override
  public String[] getParameterNames() {
    return paramNames;
  }

	
  @Override
  public Object[] getParameterValues() {
    return new Object[]{FoldX, FoldY, FoldZ, FoldW, Scale, color_mode, color_scale};
  }

    @Override
  public void setParameter(String pName, double pValue) {
    if (PARAM_FoldX.equalsIgnoreCase(pName)) 
      FoldX = pValue;
    else if (PARAM_FoldY.equalsIgnoreCase(pName)) 
      FoldY = pValue;
    else if (PARAM_FoldZ.equalsIgnoreCase(pName)) 
      FoldZ = pValue;
    else if (PARAM_FoldW.equalsIgnoreCase(pName)) 
      FoldW = pValue;
    else if (PARAM_Scale.equalsIgnoreCase(pName)) 
      Scale = pValue;
	else if (PARAM_COLOR_MODE.equalsIgnoreCase(pName))
	  color_mode = limitIntVal(Tools.FTOI(pValue), 0, 6); 
	else if (PARAM_COLOR_SCALE.equalsIgnoreCase(pName))
	  color_scale = pValue; 
	  
    else
      throw new IllegalArgumentException(pName);
  }

  @Override
  public String getName() {
    return "Naridar3DC";
  }
  
  
  }