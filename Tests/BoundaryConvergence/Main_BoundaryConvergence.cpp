#include "mpi.h"
#include "Grids.hpp"
#include "BiCGStabSolver.H"
#include "DerivativeOperators.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int result = 0;
    {
        GRParmParse pp(argc-2, argv+2, NULL, argv[1]);
        Grids::params_t p;
        Grids::read_params(pp, p);
        if (p.numLevels != 1) MayDay::Error("Uniform-grid test only");
        Grids grids(p, nullptr, false);
        grids.set_grids();
        const auto &dbl = grids.grids_data[0];
        const double dx = p.coarsestDx, L = p.domainLength[0];
        // Populate only physical cells. Every ghost starts as NaN, including
        // mixed lower-reflection/upper-extrapolation edges and corners.
        LevelData<FArrayBox> out(dbl, NUM_GRCHOMBO_VARS, 3*IntVect::Unit);
        BoundaryConditions bc;
        bc.define(dx,p.boundary_params,grids.vectDomain[0],3);
        for (DataIterator dit=dbl.dataIterator();dit.ok();++dit)
        {
            out[dit].setVal(std::numeric_limits<double>::quiet_NaN());
            for (BoxIterator bit(dbl[dit]);bit.ok();++bit)
            {
                const auto iv=bit();
                double x=(iv[0]+.5)*dx/L,y=(iv[1]+.5)*dx/L,z=(iv[2]+.5)*dx/L;
                out[dit](iv,c_lapse)=1+.01*(x*x+y*y+z*z+x*x*y*y);
                out[dit](iv,c_shift1)=x*(1+.1*x*x+.2*y*y+.3*z*z);
            }
        }
        out.exchange();
        bc.fill_grchombo_boundaries(out,Interval(c_lapse,c_shift1));
        // Check only cells outside the physical domain whose tangential
        // location belongs to this box. Intrabox ghosts are tested via exchange.
        double ghost_error=0.; int bad=0;
        for (DataIterator dit=dbl.dataIterator();dit.ok();++dit)
            for (BoxIterator bit(out[dit].box());bit.ok();++bit)
            {
                const auto iv=bit();
                if (p.coarsestDomain.domainBox().contains(iv)) continue;
                double x=(iv[0]+.5)*dx/L,y=(iv[1]+.5)*dx/L,z=(iv[2]+.5)*dx/L;
                double exact=1+.01*(x*x+y*y+z*z+x*x*y*y);
                double val=out[dit](iv,c_lapse);
                if (!std::isfinite(val)) ++bad;
                else ghost_error=std::max(ghost_error,std::abs(val-exact));
                exact=x*(1+.1*x*x+.2*y*y+.3*z*z);
                val=out[dit](iv,c_shift1);
                if (!std::isfinite(val)) ++bad;
                else ghost_error=std::max(ghost_error,std::abs(val-exact));
            }

        LevelData<FArrayBox> mg(dbl,NUM_MULTIGRID_VARS,3*IntVect::Unit);
        auto bp=p.boundary_params; bp.regularised_part_psi=1.7;
        BoundaryConditions mgbc; mgbc.define(dx,bp,grids.vectDomain[0],3);
        for (DataIterator dit=dbl.dataIterator();dit.ok();++dit)
        {
            mg[dit].setVal(std::numeric_limits<double>::quiet_NaN());
            for (BoxIterator bit(dbl[dit]);bit.ok();++bit)
                for (int c=0;c<NUM_MULTIGRID_VARS;++c)
                    mg[dit](bit(),c)=(c==c_psi_reg)?1.7:0.;
        }
        mg.exchange();
        // Both a full refresh and a solver-only refresh must preserve g=1.7.
        mgbc.fill_multigrid_boundaries(mg,Interval(0,NUM_MULTIGRID_VARS-1),false);
        mgbc.fill_multigrid_boundaries(mg,Interval(c_psi_reg,c_U_0),true);
        for (DataIterator dit=dbl.dataIterator();dit.ok();++dit)
            for (BoxIterator bit(mg[dit].box());bit.ok();++bit)
            {
                double val=mg[dit](bit(),c_psi_reg);
                if (!std::isfinite(val)) ++bad;
                else ghost_error=std::max(ghost_error,std::abs(val-1.7));
            }

        Vector<LevelData<FArrayBox>*> sol(1), rhs(1);
        sol[0]=new LevelData<FArrayBox>(dbl,NUM_CONSTRAINT_VARS,3*IntVect::Unit);
        rhs[0]=new LevelData<FArrayBox>(dbl,NUM_CONSTRAINT_VARS,IntVect::Zero);
        Vector<RefCountedPtr<LevelData<FArrayBox>>> a(1),b(1);
        a[0]=RefCountedPtr<LevelData<FArrayBox>>(new LevelData<FArrayBox>(dbl,NUM_CONSTRAINT_VARS));
        b[0]=RefCountedPtr<LevelData<FArrayBox>>(new LevelData<FArrayBox>(dbl,NUM_CONSTRAINT_VARS));
        for (DataIterator dit=dbl.dataIterator();dit.ok();++dit)
        {
            (*sol[0])[dit].setVal(0);(*rhs[0])[dit].setVal(0);
            (*a[0])[dit].setVal(0);(*b[0])[dit].setVal(1);
            for (BoxIterator bit(dbl[dit]);bit.ok();++bit)
            {
                const auto iv=bit();
                double X=1-std::pow((iv[0]+.5)*dx/L,2);
                double Y=1-std::pow((iv[1]+.5)*dx/L,2);
                double Z=1-std::pow((iv[2]+.5)*dx/L,2);
                (*rhs[0])[dit](iv,c_psi)=-.2/(L*L)*(X*Y+X*Z+Y*Z);
            }
        }
        MultilevelLinearOp<FArrayBox> op;
        op.m_num_mg_iterations=4;op.m_num_mg_smooth=4;op.m_preCondSolverDepth=-1;
        grids.define_operator(op,a,b,1.,-1.);
        BiCGStabSolver<Vector<LevelData<FArrayBox>*>> solver;
        solver.define(&op,false);solver.m_verbosity=0;solver.m_normType=2;solver.m_eps=1e-12;solver.m_imax=100;
        solver.solve(sol,rhs);
        for (DataIterator dit=dbl.dataIterator();dit.ok();++dit)
            for (BoxIterator bit(dbl[dit]);bit.ok();++bit)
                out[dit](bit(),c_lapse)=1+(*sol[0])[dit](bit(),c_psi);
        out.exchange();
        bc.fill_grchombo_boundaries(out,Interval(c_lapse,c_lapse));
        double sums[4]={0,0,0,0};
        DerivativeOperators derivs(grids.vectDx[0]);
        for (DataIterator dit=dbl.dataIterator();dit.ok();++dit)
            for (BoxIterator bit(dbl[dit]);bit.ok();++bit)
            {
                const auto iv=bit();
                double X=1-std::pow((iv[0]+.5)*dx/L,2);
                double Y=1-std::pow((iv[1]+.5)*dx/L,2);
                double Z=1-std::pow((iv[2]+.5)*dx/L,2);
                double err=out[dit](iv,c_lapse)-(1+.1*X*Y*Z);
                sums[0]+=err*err;sums[1]+=1;
                if (iv[0]==p.nCells[0]-1 || iv[1]==p.nCells[1]-1 || iv[2]==p.nCells[2]-1)
                {
                    double lap;
                    derivs.scalar_Laplacian(lap,iv,out[dit],c_lapse);
                    err=lap-(*rhs[0])[dit](iv,c_psi);
                    sums[2]+=err*err;sums[3]+=1;
                }
            }
        double total[4],gmax;int btotal;
        MPI_Allreduce(sums,total,4,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(&ghost_error,&gmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        MPI_Allreduce(&bad,&btotal,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        int rank;MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        if (rank==0) std::cout<<std::setprecision(15)<<"MANUFACTURED "<<p.nCells[0]<<" "<<std::sqrt(total[0]/total[1])<<" "<<std::sqrt(total[2]/total[3])<<" "<<gmax<<" "<<btotal<<" "<<solver.m_exitStatus<<std::endl;
        result=(btotal || gmax>1e-10 || solver.m_exitStatus!=1)?1:0;
        delete sol[0];delete rhs[0];
    }
    MPI_Finalize();return result;
}
