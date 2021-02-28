////////////////////////////////////////////////////////////////////////
// Class:       PDSPMisalignment
// Plugin Type: analyzer (art v3_05_01)
// File:        PDSPMisalignment_module.cc
//
// Generated at Sat Feb 27 12:19:53 2021 by Jacob Calcutt using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include "TVector3.h"
#include "TGraph2D.h"
#include "TH2D.h"
#include "Fit/Fitter.h"
#include "Math/Functor.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TCanvas.h"

namespace pduneana {
  class PDSPMisalignment;

  struct segparam_t {
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
    double length;

    segparam_t(double ix, double iy, double iz,
               double ivx, double ivy, double ivz,
               double ilength) :
      x(ix), y(iy), z(iz), vx(ivx), vy(ivy), vz(ivz), length(ilength) {};
    segparam_t() {};
  };

  struct distortionstruct {
    double dx;
    double dy;
    double dz;
    double roll;  // rotation about x-axis through center of APA
    double pitch; // rotation about z-axis through center of APA
    double yaw;   // rotation about y-axis through center of APA

    distortionstruct(double idx, double idy, double idz, double ir,
                     double ip, double iy) :
      dx(idx), dy(idy), dz(idz), roll(ir), pitch(ip), yaw(iy) {};
    distortionstruct() {};
  };

  struct SumDistance2 {
    // the TGraph is a data member of the object
    TGraph2D *fGraph;
  
    SumDistance2(TGraph2D *g) : fGraph(g) {}
  
    // calculate distance line-point
    double distance2(double x,double y,double z, const double *p) {
      // distance line point is D= | (xp-x0) cross  ux |
      // where ux is direction of line and x0 is a point in the line (like t = 0)
      TVector3 xp(x,y,z);
      TVector3 x0(p[0], p[2], 0. );
      TVector3 x1(p[0] + p[1], p[2] + p[3], 1. );
      TVector3 u = (x1-x0).Unit();
      double d2 = ((xp-x0).Cross(u)).Mag2();
      return d2;
    }
  
    // implementation of the function to be minimized
    double operator() (const double *par) {
      assert(fGraph != 0);
      double * x = fGraph->GetX();
      double * y = fGraph->GetY();
      double * z = fGraph->GetZ();
      int npoints = fGraph->GetN();
      double sum = 0;
      for (int i  = 0; i < npoints; ++i) {
        double d = distance2(x[i],y[i],z[i],par);
        sum += d;
      }
      //if (first) {
      //   std::cout << "Total Initial distance square = " << sum << std::endl;
      //}
      //first = false;
      return sum;
    }
  
  };


  void line(double t, const double *p, double &x, double &y, double &z) {
    // a parametric line is define from 6 parameters but 4 are independent
    // x0,y0,z0,z1,y1,z1 which are the coordinates of two points on the line
    // can choose z0 = 0 if line not parallel to x-y plane and z1 = 1;
    x = p[0] + p[1]*t;
    y = p[2] + p[3]*t;
    z = t;
  }
}


using pduneana::segparam_t;
using pduneana::distortionstruct;

class pduneana::PDSPMisalignment : public art::EDAnalyzer {
public:
  explicit PDSPMisalignment(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PDSPMisalignment(PDSPMisalignment const&) = delete;
  PDSPMisalignment(PDSPMisalignment&&) = delete;
  PDSPMisalignment& operator=(PDSPMisalignment const&) = delete;
  PDSPMisalignment& operator=(PDSPMisalignment&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  void beginJob() override;
  void endJob() override;

private:

  std::vector<TH2D *> fC2Hists;
  TGraph2D * fOutputGraph;
  TCanvas * fC1, * fC2;

  std::string fSpacePointTag;
  double fXOffset, fYOffset, fZOffset;
  double fXWidth, fYWidth, fZWidth;
  double fDistCutZ;
  double fDistZExclude;
  double fD3DClus;
  bool fPlotAllSegs, fPlotMatchedSegs;

  int fGapList1[4] = {0,1,2,3};
  int fGapList2[4] = {2,3,4,5};

  size_t GetAPAid(TVector3 & sp);
  void GetAPACenter(size_t apaid, TVector3 & apacenter);
  float DistSeg(segparam_t &seg1, segparam_t &seg2);
  bool SegMatchTest(segparam_t &seg1, segparam_t &seg2);
  void ApplyDistortion(TVector3 &sp_undistorted, distortionstruct *distlist,
                       TVector3 &sp_distorted);
  double FitSegment(std::vector<TVector3> &ptlist, std::vector<size_t>,
                    segparam_t &segparams);
  void FindMatchedSegments(std::vector<TVector3> &ptlist);
};


pduneana::PDSPMisalignment::PDSPMisalignment(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fSpacePointTag(p.get<std::string>("SpacePointTag")),
    fXOffset(p.get<double>("XOffset")),
    fYOffset(p.get<double>("YOffset")),
    fZOffset(p.get<double>("ZOffset")),
    fXWidth(p.get<double>("XWidth")),
    fYWidth(p.get<double>("YWidth")),
    fZWidth(p.get<double>("ZWidth")),
    fDistCutZ(p.get<double>("DistCutZ")),
    fDistZExclude(p.get<double>("DistZExclude")),
    fD3DClus(p.get<double>("D3DClus")),
    fPlotAllSegs(p.get<bool>("PlotAllSegs")),
    fPlotMatchedSegs(p.get<bool>("PlotMatchedSegs")) {

}

void pduneana::PDSPMisalignment::analyze(art::Event const& e) {
  std::vector<TVector3> ptlist;
  auto const& spacepoints
      = e.getValidHandle<std::vector<recob::SpacePoint>>(fSpacePointTag);

  if (!(*spacepoints).empty()) {
    for (size_t isp = 0; isp < (*spacepoints).size(); ++isp) {
        TVector3 pt( (*spacepoints)[isp].XYZ()[0],
                     (*spacepoints)[isp].XYZ()[1],
                     (*spacepoints)[isp].XYZ()[2] );
  
        // only save the points near gaps
        if (pt.Z()<200 || pt.Z()>600) continue;
  
        size_t apaid = GetAPAid(pt);
        TVector3 apacenter;
        GetAPACenter(apaid, apacenter);
        TVector3 dsp = pt - apacenter;
  
        if ( (TMath::Abs(dsp.Z()) > fZWidth/2 - fDistCutZ) &&
             (TMath::Abs(dsp.Z()) < fZWidth/2 - fDistZExclude)) 
          {
            ptlist.push_back(pt);
          }
        //std::cout << "Space point: " << x << " " << y << " " << z << endl;
    }
    // put list of space points in event selection if it passes requirements.
    FindMatchedSegments(ptlist);
  }

}

size_t pduneana::PDSPMisalignment::GetAPAid(TVector3 & sp) {
  double x = sp.X() - fXOffset;
  //double y = sp.Y() - fYOffset;
  double z = sp.Z() - fZOffset;

  size_t xindex = 0;

  // for ProtoDUNE -- we have a CPA in the middle, not an APA
  //if (x < -XWIDTH_NOMINAL) xindex = 0;
  //else if (x < XWIDTH_NOMINAL) xindex = 1;
  //else xindex = 2;
  if (x>0) xindex = 1;

  // only have one row in y

  //size_t yindex = 0;
  //if (y<0) yindex = 0;
  //else yindex = 1;

  size_t zindex = z/fZWidth;

  //std::cout << "APA id: " << xindex + 3*yindex + 6*zindex << std::endl;
  //return xindex + 3*yindex + 6*zindex;
  //std::cout << "APA id: " << xindex + 2*zindex << std::endl;
  return xindex + 2*zindex;
}

void pduneana::PDSPMisalignment::GetAPACenter(size_t apaid, TVector3 & apacenter) {
  int xindex = apaid % 2;
  //size_t yindex = ((apaid % 6) - xindex )/3;
  size_t yindex = 0;
  // was 6 for DUNE FD, now 2 for ProtoDUNE-SP
  size_t zindex = apaid / 2;

  apacenter.SetX( (2*xindex-1)*fXWidth);
  apacenter.SetY( ( ((double) yindex) - 0.5 )*fYWidth );
  apacenter.SetZ( ( ((double) zindex) + 0.5 )*fZWidth );
}

bool pduneana::PDSPMisalignment::SegMatchTest(segparam_t &seg1, segparam_t &seg2) {
  float dist1 = DistSeg(seg1,seg2);
  float dist2 = DistSeg(seg2,seg1);

  if (dist1 > 3) return false;
  if (dist2 > 3) return false;
  TVector3 p1(seg1.x,seg1.y,seg1.z);
  TVector3 p2(seg2.x,seg2.y,seg2.z);
  TVector3 v1(seg1.vx,seg1.vy,seg1.vz);
  TVector3 v2(seg2.vx,seg2.vy,seg2.vz);
  size_t apaid1 = GetAPAid(p1);
  size_t apaid2 = GetAPAid(p2);
  if (apaid1 == apaid2) return false;
  if (TMath::Abs(seg1.z - seg2.z) > 60) return false;
  if (TMath::Abs(v1.Dot(v2))<0.9) return false;

  return true;
}

float pduneana::PDSPMisalignment::DistSeg(segparam_t &seg1, segparam_t &seg2) {

  TVector3 p1(seg1.x,seg1.y,seg1.z);
  TVector3 p2(seg2.x,seg2.y,seg2.z);
  TVector3 v1(seg1.vx,seg1.vy,seg1.vz);
  TVector3 v2(seg2.vx,seg2.vy,seg2.vz);

  return ( (p2-p1).Cross(v1) ).Mag();
}

void pduneana::PDSPMisalignment::ApplyDistortion(TVector3 &sp_undistorted,
                                                 distortionstruct *distlist,
                                                 TVector3 &sp_distorted) {

  size_t apaid = GetAPAid(sp_undistorted);
  TVector3 apacenter;
  GetAPACenter(apaid, apacenter);
  distortionstruct distortion = distlist[apaid];
  TVector3 dpoffs(distortion.dx,distortion.dy,distortion.dz);

  // rotations are about the APA center.  Rotate the plane of the APA, but do
  // not distort the E field (an approximation).

  TVector3 spcenter = sp_undistorted - apacenter;
  TVector3 spcenterproj = spcenter;
  spcenterproj.SetX(0);              // use Y and Z of the point in the APA as the vector to rotate.

  spcenterproj.RotateX(distortion.roll);
  spcenterproj.RotateY(distortion.yaw);
  spcenterproj.RotateZ(distortion.pitch);

  // put back the X drift we had projected out, but add in the contributions from the rotation.

  TVector3 spcr = spcenterproj;
  spcr.SetX(spcenterproj.X() + spcenter.X());

  sp_distorted = apacenter + dpoffs + spcr;
  //std::cout << "Distortion vector: " << dpoffs.X() << " " << dpoffs.Y() << " " << dpoffs.Z() << std::endl;
  //std::cout << "Undistorted vector: " << sp_undistorted.X() << " " << sp_undistorted.Y() << " " << sp_undistorted.Z() << std::endl;
  //std::cout << "APA center: " <<  apacenter.X() << " " << apacenter.Y() << " " << apacenter.Z() << std::endl;
  //std::cout << "spcenter: " <<  spcenter.X() << " " << spcenter.Y() << " " << spcenter.Z() << std::endl;
  //std::cout << "spcenterproj: " <<  spcenterproj.X() << " " << spcenterproj.Y() << " " << spcenterproj.Z() << std::endl;
  //std::cout << "spcr: " <<  spcr.X() << " " << spcr.Y() << " " << spcr.Z() << std::endl;
  //std::cout << "Distorted vector: " << sp_distorted.X() << " " << sp_distorted.Y() << " " << sp_distorted.Z() << std::endl;
}

double pduneana::PDSPMisalignment::FitSegment(
    std::vector<TVector3> &ptlist, std::vector<size_t> cluspoint,
    segparam_t &segparams) {

  if (cluspoint.size() < 5) return -1;

  TGraph2D *gr = new TGraph2D();
  std::vector<float> zlist;

  for (size_t i=0;i<cluspoint.size();++i)
    {
      gr->SetPoint(i,ptlist[cluspoint[i]].X(),ptlist[cluspoint[i]].Y(),ptlist[cluspoint[i]].Z());
      //std::cout << "point in fit: " << ptlist[cluspoint[i]].X() << " " << ptlist[cluspoint[i]].Y() << " " << ptlist[cluspoint[i]].Z() << std::endl;
      zlist.push_back(ptlist[cluspoint[i]].Z());
    }

  ROOT::Fit::Fitter  fitter;
  
  // make the functor objet
  SumDistance2 sdist(gr);
  ROOT::Math::Functor fcn(sdist,4);
  // set the function and the initial parameter values
  double pStart[4] = {1,1,1,1};
  fitter.SetFCN(fcn,pStart);
  // set step sizes different than default ones (0.3 times parameter values)
  for (int i = 0; i < 4; ++i) fitter.Config().ParSettings(i).SetStepSize(0.01);

  bool ok = fitter.FitFCN();
  if (!ok) {
    // Error("line3Dfit","Line3D Fit failed");
    delete gr;
    return -1;
  }

  const ROOT::Fit::FitResult & result = fitter.Result();

  //std::cout << "Total final distance square " << result.MinFcnValue() << std::endl;
  //result.Print(std::cout);

  // store fit results in segmparams
  // line params x = p[0] + p[1]*t;
  //             y = p[2] + p[3]*t;
  //             z = t;

  double par0 = result.Parameter(0);
  double par1 = result.Parameter(1);
  double par2 = result.Parameter(2);
  double par3 = result.Parameter(3);

  // get min Z of the points in the fit
  std::sort(zlist.begin(),zlist.end());
  segparams.z = zlist.at(0);
  segparams.x = par0 + par1*segparams.z;
  segparams.y = par2 + par3*segparams.z;

  float dz = zlist.back() - zlist.at(0);

  TVector3 v(par1,par3,1);
  float vinorm = 1.0/v.Mag();
  v *= vinorm;
  segparams.vx = v.X();
  segparams.vy = v.Y();
  segparams.vz = v.Z();
  segparams.length = TMath::Abs(dz / segparams.vz);

  delete gr;
  return result.MinFcnValue();
}

void pduneana::PDSPMisalignment::FindMatchedSegments(
    std::vector<TVector3> &ptlist) {
  std::vector<segparam_t> seglist;
  std::vector<std::vector<size_t> > segcluslist;

  size_t numpts = ptlist.size();
  if (numpts<20) return;  // kind of minimal..

  std::vector<size_t> ptibyapa[150];   // one list of point indices per APA
  for (size_t i=0;i<numpts;++i)
    {
      size_t apaid = GetAPAid(ptlist[i]);
      ptibyapa[apaid].push_back(i);
    }

  for (size_t iapa=0;iapa<150;++iapa)
    {
      // find the clusters of points and fit lines.
      std::vector< std::vector<size_t> > cluslist;
    
      for (size_t ipoint=0;ipoint<ptibyapa[iapa].size();++ipoint)
	{
	  size_t pti = ptibyapa[iapa][ipoint];
	  // check to see if this point is close enough to any point in any existing cluster in this APA.  If so, add it.
	  bool assoc=false;
	  for (size_t iclustest = 0; iclustest<cluslist.size();++iclustest)
	    {
	      for (size_t icp=0; icp < cluslist[iclustest].size(); ++icp)
		{
		  double dist=(ptlist[pti]-ptlist[cluslist[iclustest][icp]]).Mag();
		  if (dist < fD3DClus)
		    {
		      cluslist[iclustest].push_back(pti);
		      assoc = true;
		      break;
		    }
		} // end loop over test points in cluster
	      if (assoc) break;
	    } // end loop over clusters to test
	  if (!assoc)   // if not, make a new cluster just for this point.  Todo -- merge clusters?
	    {
	      std::vector<size_t> clus;
	      clus.push_back(pti);
	      cluslist.push_back(clus);
	    }
	} // end loop over spacepoint in this APA
    
      // see if these clusters are valid segments
    
      for (size_t iclus=0; iclus<cluslist.size();++iclus)
	{
	  segparam_t segparams;
	  double fitchisq = FitSegment(ptlist,cluslist[iclus],segparams);
	  
	  bool issegment = ( fitchisq > 0  );
	  if (cluslist[iclus].size()>0)
	    { 
	      issegment &= (fitchisq/cluslist[iclus].size() < 10.0);
	    }
	  if (issegment)
	    {
	      seglist.push_back(segparams);
	      segcluslist.push_back(cluslist[iclus]);
              if (fPlotAllSegs) {
	        for (size_t i=0;i<cluslist[iclus].size();++i) {
	            int nptmp = fOutputGraph->GetN();
	            fOutputGraph->SetPoint(nptmp,ptlist[cluslist[iclus][i]].Z(),ptlist[cluslist[iclus][i]].X(),ptlist[cluslist[iclus][i]].Y());
                }
              }
	    }
	}

    } // end loop over iapa

  if (fPlotMatchedSegs) {

    for (size_t iseg=0; iseg<seglist.size()-1; ++iseg)
      {
        for (size_t jseg=iseg+1; jseg<seglist.size(); ++jseg)
          {
            if (SegMatchTest(seglist[iseg], seglist[jseg]))
              {
                for (size_t i=0;i<segcluslist[iseg].size();++i)
          	{
          	  int nptmp = fOutputGraph->GetN();
          	  fOutputGraph->SetPoint(nptmp,ptlist[segcluslist[iseg][i]].Z(),ptlist[segcluslist[iseg][i]].X(),ptlist[segcluslist[iseg][i]].Y());
          	}
                for (size_t i=0;i<segcluslist[jseg].size();++i)
          	{
          	  int nptmp = fOutputGraph->GetN();
          	  fOutputGraph->SetPoint(nptmp,ptlist[segcluslist[jseg][i]].Z(),ptlist[segcluslist[jseg][i]].X(),ptlist[segcluslist[jseg][i]].Y());
          	}

                //if (segcluslist[iseg].size() == 0) continue;
                //int iapa = GetAPAid(ptlist[segcluslist[iseg][0]]);
                //if (segcluslist[jseg].size() == 0) continue;
                //int japa = GetAPAid(ptlist[segcluslist[jseg][0]]);
                //std::cout << "APA pair: " << iapa << " " << japa << std::endl;
              }
          }
      }
  }


  // loop over parameter variations, shift space points, and fit pairs of segments, accumulating chisquareds in histograms
  // first attempt -- just try to fit for gap parameters, so loop over four gaps instead of six APA's.  When adjusting the
  // gap dimensions, just move one of the two APA's.  In general, one would have to do a simultaneous fit over all six APA's
  // positions and angles.

  for (int igap=0; igap<4; ++igap)
    {

      // fit all pairs of segments with the distortions applied
      // skip pairs of segments that aren't in the gap list
      // the joint fit will have to refit all segments however

      TH2D *htmp = (TH2D*) fC2Hists[igap]->Clone("htmp");
      htmp->Reset();
      bool hasfailed = false;

      for (size_t iseg=0; iseg<seglist.size()-1; ++iseg)
	{
	  if (segcluslist[iseg].size() == 0) continue;
	  int iapa = GetAPAid(ptlist[segcluslist[iseg][0]]);
	  if (iapa != fGapList1[igap] && iapa != fGapList2[igap]) continue;

	  for (size_t jseg=iseg+1; jseg<seglist.size(); ++jseg)
	    {
	      if (segcluslist[jseg].size() == 0) continue;
	      int japa = GetAPAid(ptlist[segcluslist[jseg][0]]);
	      //std::cout << "APA pair: " << iapa << " " << japa << std::endl;
	      if (japa != fGapList2[igap] && japa != fGapList2[igap]) continue;

	      if (SegMatchTest(seglist[iseg], seglist[jseg]) )
		{
		  for (int idx=0; idx<21; ++idx)
		    {
		      for (int idz=0; idz<21; ++ idz)
			{
			  // prepare the distortion parameters

			  double dx = -2 + ((double) idx)/5.0;
			  double dz = -2 + ((double) idz)/5.0;
			  distortionstruct distlist[6];
			  for (int iapa=0;iapa<6;++iapa)
			    {
			      distlist[iapa].dx = 0;
			      distlist[iapa].dy = 0;
			      distlist[iapa].dz = 0;
			      distlist[iapa].roll = 0;
			      distlist[iapa].pitch = 0;
			      distlist[iapa].yaw = 0;		  
			    }
			  distlist[fGapList1[igap]].dx = dx;
			  distlist[fGapList1[igap]].dz = dz;

			  std::vector<TVector3> distortedptvec;
			  // need an array of indices for all the points we want to fit -- we want to fit all of the points in the
			  // segment pair list.
			  std::vector<size_t> cpt;

			  for (size_t i=0;i<segcluslist[iseg].size();++i)
			    {
			      TVector3 ivec=ptlist[segcluslist[iseg][i]];
			      TVector3 dvec;
			      ApplyDistortion(ivec,distlist,dvec);
			      distortedptvec.push_back(dvec);
			      cpt.push_back(cpt.size());
			    }
			  for (size_t i=0;i<segcluslist[jseg].size();++i)
			    {
			      TVector3 ivec=ptlist[segcluslist[jseg][i]];
			      TVector3 dvec;
			      ApplyDistortion(ivec,distlist,dvec);
			      distortedptvec.push_back(dvec);
			      cpt.push_back(cpt.size());
			    }
			  segparam_t spar;
			  double chisquared = FitSegment(distortedptvec,cpt,spar);
			  //std::cout << "iseg, jseg, idx, idz, chisquared: " << iseg << " " << jseg << " " << idx << " " << idz << " " << chisquared << std::endl;
			  if (chisquared<0) hasfailed=true;
			  htmp->Fill(dx,dz,chisquared);
			}
		    }
		}
	    }
	}
      if (!hasfailed) fC2Hists[igap]->Add(htmp);
    }
}

void pduneana::PDSPMisalignment::beginJob() {
  gStyle->SetOptStat(0);
  gROOT->SetBatch(1);//Turns off graphics 

  art::ServiceHandle<art::TFileService> tfs;

  fOutputGraph = tfs->makeAndRegister<TGraph2D>("output_graph", "asdf");
  for (int igap=0;igap<4;++igap) {
    TString hname="igap";
    hname += igap;
    TString htitle="Chisq for gap number: ";
    htitle += igap;
    htitle += ";dx (cm); dz (cm)";
    fC2Hists.push_back(tfs->make<TH2D>(hname,htitle,21,-2.05,2.05,21,-2.05,2.05));
  }

  std::string c1name = "c";
  std::string c1title = "Event Display";
  fC1 = tfs->makeAndRegister<TCanvas>(c1name.c_str(), c1title.c_str());
//"TGraph2D Event Display",,0,0,800,800
  std::string c2name = "c2";
  std::string c2title = "Chi2squareds";
  fC2 = tfs->makeAndRegister<TCanvas>(c2name.c_str(), c2title.c_str());
//"Chisquareds",,0,0,800,800


}

void pduneana::PDSPMisalignment::endJob() {

  fOutputGraph->SetMarkerColor(1);
  fOutputGraph->SetMarkerStyle(1);
  std::string titlestring = fSpacePointTag;
  titlestring += " SpacePoints for Alignment";
  fOutputGraph->SetTitle(titlestring.c_str());
  fOutputGraph->Draw("P");
  fOutputGraph->GetXaxis()->SetTitle("Z (cm)");
  fOutputGraph->GetXaxis()->SetTitleColor(4);
  fOutputGraph->GetYaxis()->SetTitle("X (cm)");
  fOutputGraph->GetYaxis()->SetTitleColor(4);
  fOutputGraph->GetZaxis()->SetTitle("Y (cm)");
  fOutputGraph->GetZaxis()->SetTitleColor(4);
  fC1->Update();  
  //c->Print("space_points.pdf");

  fC2->Divide(2,2);

  for (int igap=0; igap<4; ++igap)
    {
      fC2->cd(igap+1);
      fC2Hists[igap]->Draw("colz");
      fC2Hists[igap]->Write();

      int iminx=0;
      int iminy=0;
      int iminz=0;
      fC2Hists[igap]->GetMinimumBin(iminx,iminy,iminz);
      std::cout << "dx, dz fit for gap: " << igap << " " << fC2Hists[igap]->GetXaxis()->GetBinCenter(iminx) << " " << fC2Hists[igap]->GetYaxis()->GetBinCenter(iminy) << std::endl;

    }
}

DEFINE_ART_MODULE(pduneana::PDSPMisalignment)
