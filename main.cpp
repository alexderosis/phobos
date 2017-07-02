///------------------------------------------------------------------------------------------------------------------------------------------------------
// This program is part of the paper "Non-orthogonal central moments 
// relaxing to a discrete equilibrium: a D2Q9 lattice Boltzmann model".
// It allows us to run the test involving the lid-driven cavity.
// ---> Please modify "nx" and "ny" in order to set different grid size.
// ---> The quantity "u_lid" is the velocity of the lid.
// ---> Non-unitary relaxation frequencies vary according to the desired
// ---> value of "Reynolds".
// ---> The number of desired time steps is "nsteps".
// ---> The program produces a .vtk file every "n_out" time steps.
// Copyright (C) 2016  Alessandro De Rosis (alessandro@bm.technion.ac.il)
// 
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
// 
// This is the "LidDrivenCavity_CentralMoments.cpp" file.
// Author: Alessandro De Rosis (alessandro@bm.technion.ac.il)
// Address: Biofluids Laboratory, Department of Biomedical Engineering, 
// Julius Silver Building, Office 254,
// Technion - Israel Institute of Technology, 32000 Haifa, Israel.
// Day: 1st November 2016
// Website: www.felbaproject.com
///------------------------------------------------------------------------------------------------------------------------------------------------------
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cassert>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;
///------------------------------------------------------------------------------------------------------------------------------------------------------
const int nx = 201, ny = 401, np = 9, nsteps = 1000000001, n_out = 2500, cx[np] = {0., 1., 0., -1., 0., 1., -1., -1., 1.}, cy[np] = {0., 0., 1., 0., -1., 1., 1., -1., -1.};
const double w[np] = {4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.}, cs2 = 1./3., rho0 = 1., u_lid = 0.15, Reynolds = 1./0., ni = 0., tau = ni*3.+.5, omega = 1./tau, omega1 = 1.-omega;
double f1[np][nx][ny], f2[np][nx][ny], u[nx][ny], v[nx][ny], rho[nx][ny], collision_operator[np], feq[np], temp_pop[np], k3, k4, k5, k6, k7, k8, k3_eq, k4_eq, k5_eq, k6_eq, k7_eq, k8_eq, CX, CY, U, V, R, ftemp, U2, V2, A, B, C, Pxx, Pxy, Pyy, fneq, QPI, UV, U2V2;
int newi, newj, i, j, k;
const bool plot_vtk = true;
///------------------------------------------------------------------------------------------------------------------------------------------------------
///------------------------------------------------------------------------------------------------------------------------------------------------------
void write_fluid_vtk(int time)
{
	stringstream output_filename;
	output_filename << "vtk_fluid/fluid_t" << time << ".vtk";
	ofstream output_file;

	output_file.open(output_filename.str().c_str());

	output_file << "# vtk DataFile Version 3.0\n";
	output_file << "fluid_state\n";
	output_file << "ASCII\n";
	output_file << "DATASET RECTILINEAR_GRID\n";
	output_file << "DIMENSIONS " << nx << " " << ny << " 1" << "\n";
	output_file << "X_COORDINATES " << nx << " float\n";
	for(int i = 0; i < nx; ++i)
		output_file << i << " ";
	output_file << "\n";
	output_file << "Y_COORDINATES " << ny  << " float\n";
	for(int j = 0; j < ny ; ++j)
		output_file << j  << " ";
	output_file << "\n";
	output_file << "Z_COORDINATES " << 1 << " float\n";
	output_file << 0 << "\n";
	output_file << "POINT_DATA " << (nx) * (ny) << "\n";

	output_file << "VECTORS velocity_vector float\n";
	for(int Y = 0; Y < ny ; ++Y)
		for(int X = 0; X < nx; ++X)
			output_file << u[X][Y]/u_lid << " " << v[X][Y]/u_lid << " 0\n";

	output_file.close();
}
///------------------------------------------------------------------------------------------------------------------------------------------------------
void initial_state()
{
	for(i=0; i<nx; i++)
        for(j=0; j<ny; j++)
		{
			R = rho[i][j] = rho0;
			U = u[i][j] = 0;
            V = v[i][j] = 0.;
			C = -1.5*(U*U+V*V);
			for(k=0; k<np;k++)
			{
           		A = U*cx[k]+V*cy[k];
           		B = 4.5*A*A;
                f1[k][i][j] = f2[k][i][j] = w[k]*R*(1.+3.*A+B+C);
			}
		}
}
///------------------------------------------------------------------------------------------------------------------------------------------------------
int my_algorithm()
{
	int check = 0;
	for(i=0; i<nx; i++)
		for(j=0; j<ny; j++)
		{
			/// compute macroscopic variables
			R = U = V = 0.;
			for(k=0; k<np; k++)
			{
				ftemp = temp_pop[k] = f1[k][i][j];
				R += ftemp;
				U += ftemp*cx[k];
				V += ftemp*cy[k];
			}
			U /= R;
			V /= R;
			u[i][j] = U;
			v[i][j] = V;
			rho[i][j] = R;
			if(fabs(U)>1)
				check = 1;
			U2 = U*U;
			V2 = V*V;
			UV = U*V;
			U2V2 = U2*V2;
			/// compute pre-collision central moments
			k4 = k5 = 0.;
			for(k=0; k<np; k++)
			{	
				CX = cx[k]-U;
				CY = cy[k]-V;			
				ftemp = temp_pop[k];
				// k3 += ftemp*CX*CX+ftemp*CY*CY;
				k4 += ftemp*(CX*CX-CY*CY);
				k5 += ftemp*CX*CY;
				// k6 += ftemp*CX*CX*CY;
				// k7 += ftemp*CX*CY*CY;
				// k8 += ftemp*CX*CX*CY*CY;
			}
			/// evaluate equilibrium central moments (see CentralMoments_D2Q9.m)
			k3_eq = R*2*cs2;
            k6_eq = -R*U2*V;
            k7_eq = -R*U*V2;
 			k8_eq = R*(1+27*U2V2)*cs2*cs2;
 			/// post-collision state
			k3 = k3_eq;
			k4 = omega1*k4;
			k5 = omega1*k5;
			k6 = k6_eq;
			k7 = k7_eq;				
			k8 = k8_eq;
			/// reconstruct post-collision populations
			collision_operator[0] = k3*(-1+0.5*U2+0.5*V2)+k4*(-U2+V2)+4*UV*k5+2*V*k6+2*U*k7+k8+R*(1-U2-V2+U2V2);
    		collision_operator[1] = 0.25*k3*(1-U-U2-V2)+0.25*k4*(1+U+U2-V2)+k5*(-V-2*UV)-V*k6+k7*(-0.5-U)-0.5*k8+0.5*R*(U+U2-U*V2-U2V2);
    		collision_operator[2] = 0.25*k3*(1-V-U2-V2)+0.25*k4*(-1-V+U2-V2)+k5*(-U-2*UV)+k6*(-0.5-V)-U*k7-0.5*k8+0.5*R*(V+V2-U2*V-U2V2);
	    	collision_operator[3] = 0.25*k3*(1+U-U2-V2)+0.25*k4*(1-U+U2-V2)+k5*(V-2*UV)-V*k6+k7*(0.5-U)-0.5*k8+0.5*R*(-U+U2+U*V2-U2V2);
	    	collision_operator[4] = 0.25*k3*(1+V-U2-V2)+0.25*k4*(-1+V+U2-V2)+k5*(U-2*UV)+k6*(0.5-V)-U*k7-0.5*k8+0.5*R*(-V+V2+U2*V-U2V2);
	    	collision_operator[5] = 0.125*k3*(U+V+U2+V2)+0.125*k4*(-U+V-U2+V2)+k5*(0.25+0.5*U+0.5*V+UV)+k6*(0.25+0.5*V)+k7*(0.25+0.5*U)+0.25*k8+0.25*R*(UV+U2*V+U*V2+U2V2);
	    	collision_operator[6] = 0.125*k3*(-U+V+U2+V2)+0.125*k4*(U+V-U2+V2)+k5*(-0.25+0.5*U-0.5*V+UV)+k6*(0.25+0.5*V)+k7*(-0.25+0.5*U)+0.25*k8+0.25*R*(-UV+U2*V-U*V2+U2V2);
	    	collision_operator[7] = 0.125*k3*(-U-V+U2+V2)+0.125*k4*(U-V-U2+V2)+k5*(0.25-0.5*U-0.5*V+UV)+k6*(-0.25+0.5*V)+k7*(-0.25+0.5*U)+0.25*k8+0.25*R*(UV-U2*V-U*V2+U2V2);
	    	collision_operator[8] = 0.125*k3*(U-V+U2+V2)+0.125*k4*(-U-V-U2+V2)+k5*(-0.25-0.5*U+0.5*V+UV)+k6*(-0.25+0.5*V)+k7*(0.25+0.5*U)+0.25*k8+0.25*R*(-UV-U2*V+U*V2+U2V2);
			/// save the post-collision state in f1 and stream in f2
			for(k=0; k<np; k++)
			{
				f1[k][i][j] = collision_operator[k];
				newi = i+cx[k];
				newj = j+cy[k];
				if(i==0 || i==nx-1)
					newi = (newi+nx)%nx;
				if(j==0 || j==ny-1)
					newj = (newj+ny)%ny;
				f2[k][newi][newj] = collision_operator[k];
			}

		}
	return check;
}
///------------------------------------------------------------------------------------------------------------------------------------------------------
void boundary_conditions() /// by regularized scheme
{
	/// west and east walls
	for(j=0; j<ny; j++)
	{
    	Pxx = Pxy = Pyy = 0.;
	    V = v[0][j] = 0.;
       	U = u[0][j] = 0;
       	R = rho[0][j] = rho0;
       	C = -1.5*(U*U+V*V);
       	for(k=0; k<np; k++)
       	{
           	A = U*cx[k]+V*cy[k];
           	B = 4.5*A*A;
           	feq[k] = w[k]*R*(1.+3.*A+B+C);
           	fneq = f1[k][0][j]-feq[k];
           	Pxx += fneq*cx[k]*cx[k];
           	Pxy += fneq*cx[k]*cy[k];
           	Pyy += fneq*cy[k]*cy[k];
       	}
       	for(k=0; k<np; k++)
       	{		
           	QPI = Pxx*(cx[k]*cx[k]-cs2)+2*Pxy*cx[k]*cy[k]+Pyy*(cy[k]*cy[k]-cs2);
           	f2[k][0][j] = feq[k]+4.5*w[k]*QPI;
       	}
       	Pxx = Pxy = Pyy = 0.;
	    V = v[nx-1][j] = 0.;
       	U = u[nx-1][j] = 0.;
       	R = rho[nx-1][j] = rho0;
       	C = -1.5*(U*U+V*V);
       	for(k=0; k<np; k++)
       	{
           	A = U*cx[k]+V*cy[k];
           	B = 4.5*A*A;
           	feq[k] = w[k]*R*(1.+3.*A+B+C);
           	fneq = f1[k][nx-1][j]-feq[k];
           	Pxx += fneq*cx[k]*cx[k];
           	Pxy += fneq*cx[k]*cy[k];
           	Pyy += fneq*cy[k]*cy[k];
       	}
       	for(k=0; k<np; k++)
       	{
           	QPI = Pxx*(cx[k]*cx[k]-cs2)+2*Pxy*cx[k]*cy[k]+Pyy*(cy[k]*cy[k]-cs2);
           	f2[k][nx-1][j] = feq[k]+4.5*w[k]*QPI;
       	}
	}
	/// south wall and north lid
	for(i=0; i<nx; i++)
	{
   		Pxx = Pxy = Pyy = 0.;
	    V = v[i][0] = 0.;
       	U = u[i][0] = 0.;
       	R = rho[i][0] = rho0;
       	C = -1.5*(U*U+V*V);
       	for(k=0; k<np; k++)
       	{
           	A = U*cx[k]+V*cy[k];
           	B = 4.5*A*A;
           	feq[k] = w[k]*R*(1.+3.*A+B+C);
           	fneq = f1[k][i][0]-feq[k];
           	Pxx += fneq*cx[k]*cx[k];
           	Pxy += fneq*cx[k]*cy[k];
           	Pyy += fneq*cy[k]*cy[k];
       	}
       	for(k=0; k<np; k++)
       	{		
           	QPI = Pxx*(cx[k]*cx[k]-cs2)+2*Pxy*cx[k]*cy[k]+Pyy*(cy[k]*cy[k]-cs2);
           	f2[k][i][0] = feq[k]+4.5*w[k]*QPI;
       	}
       	Pxx = Pxy = Pyy = 0.;
	    V = v[i][ny-1] = 0.;
       	U = u[i][ny-1] = u_lid;
       	R = rho[i][ny-1] = rho0;
       	C = -1.5*(U*U+V*V);
       	for(k=0; k<np; k++)
       	{
           	A = U*cx[k]+V*cy[k];
           	B = 4.5*A*A;
           	feq[k] = w[k]*R*(1.+3.*A+B+C);
           	fneq = f1[k][i][ny-1]-feq[k];
           	Pxx += fneq*cx[k]*cx[k];
           	Pxy += fneq*cx[k]*cy[k];
           	Pyy += fneq*cy[k]*cy[k];
       	}
       	for(k=0; k<np; k++)
       	{
           	QPI = Pxx*(cx[k]*cx[k]-cs2)+2*Pxy*cx[k]*cy[k]+Pyy*(cy[k]*cy[k]-cs2);
           	f2[k][i][ny-1] = feq[k]+4.5*w[k]*QPI;
       	}
	 }
}
///------------------------------------------------------------------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
	int check_convergence, t;
	system("mkdir vtk_fluid");
	initial_state();
	for(t=0; t<nsteps; t++)
    {
		check_convergence = my_algorithm(); // compute macroscopic variables, collide and stream
		boundary_conditions(); // apply boundary conditions 
		memcpy(f1, f2, sizeof(f2)); // swap f2 in f1
		if(plot_vtk==true && t%n_out==0)
			write_fluid_vtk(t);
		if(check_convergence==1 && t>1)
			goto labelA;
    }
    labelA:
    return 0;
}
///------------------------------------------------------------------------------------------------------------------------------------------------------
///------------------------------------------------------------------------------------------------------------------------------------------------------
