#include <iostream>
#include <cucheb.h>
using namespace std;

#ifndef hamop_h_included
#define hamop_h_included

class Hamop{
	private:
		int dims;
		int nx;
		int ny;
		int nz;
		double xmin, xmax;
		double ymin, ymax;
		double zmin, zmax;
		double xpot, ypot, zpot;

	public:

		// sparse matrix routines
		void Dcsr(double *vals, int *rows, int *cols, double shift);
		void Dcsc(double *vals, int *rows, int *cols, double shift);
		void DcscLowerTri(double *vals, int *rows, int *cols, double shift);

		// matrix vector product
		void Mult(double* vecin, double* vecout);

		// constructors
		Hamop(void):dims(1),nx(1),ny(1),nz(1),xmin(-1.0),xmax(1.0),ymin(-1.0),ymax(1.0),zmin(-1.0),zmax(1.0),xpot(1.0),
			ypot(1.0),zpot(1.0){}
		Hamop(const Hamop * ho){
			setDims(ho->getDims());
			setNx(ho->getNx());
			setNy(ho->getNy());
			setNz(ho->getNz());
			setXminXmax(ho->getXmin(),ho->getXmax());
			setYminYmax(ho->getYmin(),ho->getYmax());
			setZminZmax(ho->getZmin(),ho->getZmax());
			setXpotYpotZpot(ho->getXpot(),ho->getYpot(),ho->getZpot());
		}

		// printer
		void print(void){
			std::cout<< "Hamop:" << "\n";
			std::cout<< " dims = " << dims << "\n";
			std::cout<< " nx = " << nx << "\n";
			std::cout<< " ny = " << ny << "\n";
			std::cout<< " nz = " << nz << "\n";
			std::cout<< " [xmin,xmax] = [" << xmin << "," << xmax << "]\n";
			std::cout<< " [ymin,ymax] = [" << ymin << "," << ymax << "]\n";
			std::cout<< " [zmin,zmax] = [" << zmin << "," << zmax << "]\n";
			std::cout<< " (xpot,ypot,zpot) = (" << xpot << "," << ypot << "," << zpot << ")\n";
			std::cout<< "\n";
		}

		// gets
		int getDims(void) const{return dims;}
		int getNx(void) const{return nx;}
		int getNy(void) const{return ny;}
		int getNz(void) const{return nz;}
		double getXmin(void) const{return xmin;}
		double getXmax(void) const{return xmax;}
		double getYmin(void) const{return ymin;}
		double getYmax(void) const{return ymax;}
		double getZmin(void) const{return zmin;}
		double getZmax(void) const{return zmax;}
		double getXpot(void) const{return xpot;}
		double getYpot(void) const{return ypot;}
		double getZpot(void) const{return zpot;}

		// sets
		void setDims(int D){
			if(D == 1 || D == 2 || D == 3){
				dims = D;
			}
			else{
				dims = 1;
				std::cout<< "\nWarning in Hamop:" << "\n";
				std::cout<< " Not a valid setting for dims!" << "\n";
				std::cout<< " Options: 1, 2, 3" << "\n";
				std::cout<< " Setting dims to 1." << "\n\n";
			}
		}
		void setNx(int N){
			if(N > 0){
				nx = N;
			}
			else{
				nx = 1;
				std::cout<< "\nWarning in Hamop:" << "\n";
				std::cout<< " Not a valid setting for nx!" << "\n";
				std::cout<< " Must be a positive integer." << "\n";
				std::cout<< " Setting nx to 1." << "\n\n";
			}
		}
		void setNy(int N){
			if(N > 0){
				ny = N;
			}
			else{
				ny = 1;
				std::cout<< "\nWarning in Hamop:" << "\n";
				std::cout<< " Not a valid setting for ny!" << "\n";
				std::cout<< " Must be a positive integer." << "\n";
				std::cout<< " Setting ny to 1." << "\n\n";
			}
		}
		void setNz(int N){
			if(N > 0){
				nz = N;
			}
			else{
				nz = 1;
				std::cout<< "\nWarning in Hamop:" << "\n";
				std::cout<< " Not a valid setting for nz!" << "\n";
				std::cout<< " Must be a positive integer." << "\n";
				std::cout<< " Setting nz to 1." << "\n\n";
			}
		}
		void setXminXmax(double A, double B){
			if(B > A){
				xmin = A;
				xmax = B;
			}
			else{
				xmin = -1.0;
				xmax = 1.0;
				std::cout<< "\nWarning in Hamop:" << "\n";
				std::cout<< " Not a valid setting for xmin and xmax!" << "\n";
				std::cout<< " Must be a real interval with non-empty interior." << "\n";
				std::cout<< " Setting [xmin,xmax] to [-1,1]." << "\n\n";
			}
		}
		void setYminYmax(double A, double B){
			if(B > A){
				ymin = A;
				ymax = B;
			}
			else{
				ymin = -1.0;
				ymax = 1.0;
				std::cout<< "\nWarning in Hamop:" << "\n";
				std::cout<< " Not a valid setting for ymin and ymax!" << "\n";
				std::cout<< " Must be a real interval with non-empty interior." << "\n";
				std::cout<< " Setting [ymin,ymax] to [-1,1]." << "\n\n";
			}
		}
		void setZminZmax(double A, double B){
			if(B > A){
				zmin = A;
				zmax = B;
			}
			else{
				zmin = -1.0;
				zmax = 1.0;
				std::cout<< "\nWarning in Hamop:" << "\n";
				std::cout<< " Not a valid setting for zmin and zmax!" << "\n";
				std::cout<< " Must be a real interval with non-empty interior." << "\n";
				std::cout<< " Setting [zmin,zmax] to [-1,1]." << "\n\n";
			}
		}
		void setXpotYpotZpot(double A, double B, double C){
			if(A >= 0 && B >= 0 && C >= 0){
				xpot = A;
				ypot = B;
				zpot = C;
			}
			else{
				xpot = 1.0;
				ypot = 1.0;
				zpot = 1.0;
				std::cout<< "\nWarning in Hamop:" << "\n";
				std::cout<< " Not a valid setting for xpot, ypot, and zpot!" << "\n";
				std::cout<< " Must be non-negative, real numbers." << "\n";
				std::cout<< " Setting xpot = ypot = zpot = 1." << "\n\n";
			}
		}
};


// double precision subroutines
void DmultHamop(double* vecin, double* vecout);
void DmultLap(int dim, int nx, int ny, int nz, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, const double* vecin, double* vecout,int add);
__global__ void Dlap(int dim, int nx, int ny, int nz, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, const double* vecin, double* vecout,int add);
void DmultPot(int dims, int nx, int ny, int nz, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double xpot, double ypot, double zpot, const double* vecin, double* vecout,int add);
__global__ void Dpot(int dims, int nx, int ny, int nz, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double xpot, double ypot, double zpot, const double* vecin, double* vecout,int add);

#endif
