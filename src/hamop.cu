#include <hamop.h>

// Double Precision GPU Finite Difference Hamiltonian subroutines
void Hamop::Mult(double* vecin, double* vecout){
	
	// compute variables
	int n = nx*ny*nz;
	double scl = 1.0;

	// along x dim
	DmultLap(1,nx,ny,nz,xmin,xmax,ymin,ymax,zmin,zmax,vecin,vecout,0);

	// along y dim
	if(dims > 1){
		DmultLap(2,nx,ny,nz,xmin,xmax,ymin,ymax,zmin,zmax,vecin,vecout,1);
	}

	// along z dim
	if(dims > 2){
		DmultLap(3,nx,ny,nz,xmin,xmax,ymin,ymax,zmin,zmax,vecin,vecout,1);
	}

	// x^2 potential
	DmultPot(dims,nx,ny,nz,xmin,xmax,ymin,ymax,zmin,zmax,xpot,ypot,zpot,vecin,vecout,1);
}
void DmultLap(int dim, int nx, int ny, int nz, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, const double* vecin, double* vecout, int add){

	// compute variables
	int BLOCK_SIZE, GRID_SIZE;
	int n = nx*ny*nz;

	// set BLOCK_SIZE and compute GRID_SIZE
	BLOCK_SIZE = 512;
	if(n%BLOCK_SIZE == 0){
		GRID_SIZE = n/512;
	}
	else{
		GRID_SIZE = n/512 + 1;
	}

	// call kernel
	Dlap<<<GRID_SIZE,BLOCK_SIZE>>>(dim,nx,ny,nz,xmin,xmax,ymin,ymax,zmin,zmax,vecin,vecout,add);

}

__global__ void Dlap(int dim, int nx, int ny, int nz, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, const double* vecin, double* vecout, int add){
	
	// compute variables
	int tix = threadIdx.x, bix = blockIdx.x, bdx = blockDim.x;
	int ii = bix*bdx+tix;
	int xp, yp, zp, n;
	double scl = 2.0;
	double hx, hy, hz;

	// set variables
	n = nx*ny*nz;
	hx = (double)(nx+1)/(xmax - xmin);
	hy = (double)(ny+1)/(ymax - ymin);
	hz = (double)(nz+1)/(zmax - zmin);

	if(ii < n){
		// get pos
		xp = ii%nx;
		yp = ((ii-xp)/nx)%ny;
		zp = (((ii-xp)/nx-yp)/ny);

		if(add == 0){
			// along x dim
			if(dim == 1){
				if(xp%nx == 0){vecout[ii] = hx*hx*(scl*vecin[ii] - vecin[ii+1]);}
				else if(xp%nx == nx-1){vecout[ii] = hx*hx*(scl*vecin[ii] - vecin[ii-1]);}
				else{vecout[ii] = hx*hx*(-vecin[ii-1] + scl*vecin[ii] - vecin[ii+1]);}
			}

			// along y dim
			if(dim == 2){
				if(yp%ny == 0){vecout[ii] = hy*hy*(scl*vecin[ii] - vecin[ii+nx]);}
				else if(yp%ny == ny-1){vecout[ii] = hy*hy*(scl*vecin[ii] - vecin[ii-nx]);}
				else{vecout[ii] = hy*hy*(-vecin[ii-nx] + scl*vecin[ii] - vecin[ii+nx]);}
			}

			// along z dim
			if(dim == 3){
				if(zp%nz == 0){vecout[ii] = hz*hz*(scl*vecin[ii] - vecin[ii+nx*ny]);}
				else if(zp%nz == nz-1){vecout[ii] = hz*hz*(scl*vecin[ii] - vecin[ii-nx*ny]);}
				else{vecout[ii] = hz*hz*(-vecin[ii-nx*ny] + scl*vecin[ii] - vecin[ii+nx*ny]);}
			}
		}
		
		else{
			// along x dim
			if(dim == 1){
				if(xp%nx == 0){vecout[ii] += hx*hx*(scl*vecin[ii] - vecin[ii+1]);}
				else if(xp%nx == nx-1){vecout[ii] += hx*hx*(scl*vecin[ii] - vecin[ii-1]);}
				else{vecout[ii] += hx*hx*(-vecin[ii-1] + scl*vecin[ii] - vecin[ii+1]);}
			}

			// along y dim
			if(dim == 2){
				if(yp%ny == 0){vecout[ii] += hy*hy*(scl*vecin[ii] - vecin[ii+nx]);}
				else if(yp%ny == ny-1){vecout[ii] += hy*hy*(scl*vecin[ii] - vecin[ii-nx]);}
				else{vecout[ii] += hy*hy*(-vecin[ii-nx] + scl*vecin[ii] - vecin[ii+nx]);}
			}

			// along z dim
			if(dim == 3){
				if(zp%nz == 0){vecout[ii] += hz*hz*(scl*vecin[ii] - vecin[ii+nx*ny]);}
				else if(zp%nz == nz-1){vecout[ii] += hz*hz*(scl*vecin[ii] - vecin[ii-nx*ny]);}
				else{vecout[ii] += hz*hz*(-vecin[ii-nx*ny] + scl*vecin[ii] - vecin[ii+nx*ny]);}
			}
		}
	}
}

void DmultPot(int dims, int nx, int ny, int nz, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double xpot, double ypot, double zpot, const double* vecin, double* vecout, int add){

	// compute variables
	int BLOCK_SIZE, GRID_SIZE;
	int n = nx*ny*nz;

	// set BLOCK_SIZE and compute GRID_SIZE
	BLOCK_SIZE = 512;
	if(n%BLOCK_SIZE == 0){
		GRID_SIZE = n/512;
	}
	else{
		GRID_SIZE = n/512 + 1;
	}

	// call kernel
	Dpot<<<GRID_SIZE,BLOCK_SIZE>>>(dims,nx,ny,nz,xmin,xmax,ymin,ymax,zmin,zmax,xpot,ypot,zpot,vecin,vecout,add);

}

__global__ void Dpot(int dims, int nx, int ny, int nz, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double xpot, double ypot, double zpot, const double* vecin, double* vecout, int add){
	
	// compute variables
	int tix = threadIdx.x, bix = blockIdx.x, bdx = blockDim.x;
	int ii = bix*bdx+tix;
	int xp, yp, zp, n;
	double hx, hy, hz;
	double x, y, z;

	// set variables
	n = nx*ny*nz;
	hx = (double)(nx+1)/(xmax - xmin);
	hy = (double)(ny+1)/(ymax - ymin);
	hz = (double)(nz+1)/(zmax - zmin);

	// get pos
	xp = ii%nx;
	yp = ((ii-xp)/nx)%ny;
	zp = (((ii-xp)/nx-yp)/ny);

	// set variables
	x = xmin + (double)(xp+1)/hx;
	y = ymin + (double)(yp+1)/hy;
	z = zmin + (double)(zp+1)/hz;

	if(add == 0){
		// 1 dim
		if(dims == 1 && ii < n){vecout[ii] = (xpot*xpot*x*x)*vecin[ii];}

		// 2 dim
		if(dims == 2 && ii < n){vecout[ii] = (xpot*xpot*x*x + ypot*ypot*y*y)*vecin[ii];}

		// 3 dim
		if(dims == 3 && ii < n){vecout[ii] = (xpot*xpot*x*x + ypot*ypot*y*y + zpot*zpot*z*z)*vecin[ii];}
	}
	else{
		// 1 dim
		if(dims == 1 && ii < n){vecout[ii] += (xpot*xpot*x*x)*vecin[ii];}

		// 2 dim
		if(dims == 2 && ii < n){vecout[ii] += (xpot*xpot*x*x + ypot*ypot*y*y)*vecin[ii];}

		// 3 dim
		if(dims == 3 && ii < n){vecout[ii] += (xpot*xpot*x*x + ypot*ypot*y*y + zpot*zpot*z*z)*vecin[ii];}	
	}
}

