///////////////////////////////////////////////////////////////////////////////
//
// Author: JStock
// \brief: Stealing methods from Aleena's Michel Analysis module for use
// in my own work with stop muons.
//
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
////////////////////Stop Muon Library Headers and Classes//////////////////////
///////////////////////////////////////////////////////////////////////////////

#include "MuonFinderLib.h"


namespace MFL{
  protoana::ProtoDUNETrackUtils trackUtil;
}

bool MFL::IsPointInBounds(TVector3 const &P)
{
  auto const & v = mfl_conf.activeBounds_eff;
  auto const & fTSV = mfl_conf.fThicknessStartVolume;
  //Checks from Aleena's module.
  bool c1 = (P.Y()>=(v[3]-fTSV) && P.Y()<=v[3]);
  bool c2 = (P.X()>=v[0] && P.X()<=(v[0]+fTSV));
  bool c3 = (P.X()<=v[1] && P.X()>=(v[1]-fTSV));
  bool c4 = (P.Z()>=v[4] && P.Z()<=(v[4]+fTSV));
  bool c5 = (P.Z()<=v[5] && P.Z()>=(v[5]-fTSV));

  return c1 || c2 || c3 || c4 || c5 ;
}

bool MFL::IsPointInFV(TVector3 const & t)
{
  auto const & v = mfl_conf.fiducialBounds;
  bool c1 = t.X() >= v[0];
  bool c2 = t.X() <= v[1];
  bool c3 = t.Y() >= v[2];
  bool c4 = t.Y() <= v[3];
  bool c5 = t.Z() >= v[4];
  bool c6 = t.Z() <= v[5];
  return c1 && c2 && c3 && c4 && c5 && c6;
}

void MFL::SwitchDirection(TVector3 &ps_s, TVector3 &ps_e, TVector3 &dr_s, TVector3 dr_e)
{
  double x1(ps_s.X()),x2(ps_e.X()),y1(ps_s.Y()),y2(ps_e.Y()),z1(ps_s.Z()),z2(ps_e.Z());
  double dx1(dr_s.X()),dx2(dr_e.X()),dy1(dr_s.Y()),dy2(dr_e.Y()),dz1(dr_s.Z()),dz2(dr_e.Z());
  ps_s.SetXYZ(x2,y2,z2);
  ps_e.SetXYZ(x1,y1,z1);
  dr_s.SetXYZ(dx2,dy2,dz2);
  dr_e.SetXYZ(dx1,dy1,dz1);
}

void MFL::ForceDirection(TVector3 &ps_s, TVector3 &ps_e, TVector3 &dr_s, TVector3 dr_e)
{
  /*Check if switch needed*/
  if(MFL::IsPointInBounds(ps_e))
  {
    SwitchDirection(ps_s, ps_e, dr_s, dr_e);
  }
}

bool MFL::CathCross(TVector3 const& v1, TVector3 const& v2)
{
  return ( (v1.X()*v2.X())<0);
}

bool MFL::FidCut(TVector3 const &p)
{
  bool c1 = IsPointInFV(p);
  bool c2 = p.Y()>mfl_conf.y_end_cut;
  bool c3 = p.Z()>mfl_conf.z_end_min_cut;
  bool c4 = p.Z()<mfl_conf.z_end_max_cut;
  return c1&&c2&&c3&&c4;
}

bool MFL::BrokenTrack(art::Ptr<recob::Track> const& ptrack, std::vector<art::Ptr<recob::Track>> const& tracklist )
{
  for(size_t i=0; i<tracklist.size(); ++i)
  {
    art::Ptr<recob::Track> ptrack2 = tracklist.at(i);
    if (ptrack==ptrack2) continue;
    TVector3 ps_s(
        ptrack->LocationAtPoint(ptrack->FirstValidPoint()).X(),
        ptrack->LocationAtPoint(ptrack->FirstValidPoint()).Y(),
        ptrack->LocationAtPoint(ptrack->FirstValidPoint()).Z()
        ); //inherited from Aleena
    TVector3 ps_e = ptrack->End<TVector3>();
    TVector3 dr_s = ptrack->DirectionAtPoint<TVector3>(ptrack->FirstValidPoint());
    TVector3 dr_e = ptrack->DirectionAtPoint<TVector3>(ptrack->LastValidPoint());

    TVector3 ps_s2(
        ptrack2->LocationAtPoint(ptrack2->FirstValidPoint()).X(),
        ptrack2->LocationAtPoint(ptrack2->FirstValidPoint()).Y(),
        ptrack2->LocationAtPoint(ptrack2->FirstValidPoint()).Z()
        );
    TVector3 ps_e2 = ptrack2->End<TVector3>();
    TVector3 dr_s2(ptrack2->DirectionAtPoint<TVector3>(ptrack2->FirstValidPoint()));
    TVector3 dr_e2(ptrack2->DirectionAtPoint<TVector3>(ptrack2->LastValidPoint()));

    double bwdiststart = sqrt(
        pow(ps_s2.Y()-ps_e.Y(),2) +
        pow(ps_s2.Z()-ps_e.Z(),2)); //since the non-T0 tagged tracks have wrong x position associated with them
    double bwdistend   = sqrt(
        pow(ps_e2.Y() - ps_e.Y(),2) +
        pow(ps_e2.Z() - ps_e.Z(),2) );

    double cosopeningangleStart =
      ( cos( dr_e.Theta() ) * cos( dr_s2.Theta() ) )
      +
      ( sin(dr_e.Theta()) * sin(dr_e2.Theta()) * cos(dr_e2.Phi()-dr_e.Phi()));
    double cosopeningangleEnd   =
      ( cos( dr_e.Theta() ) * cos( dr_e2.Theta() ) )
      +
      ( sin( dr_e.Theta() ) * sin( dr_e2.Theta() ) * cos( dr_e2.Phi() - dr_e.Phi() ) ) ;
    if(
        ( abs(bwdiststart)<30 && abs(cosopeningangleStart) > 0.97)
        ||
        (abs(bwdistend)<30 && abs(cosopeningangleEnd)>0.97)
        ||
        (abs(bwdiststart)<50 && abs(cosopeningangleStart)>0.998)
        ||
        (abs(bwdistend)<50 && abs(cosopeningangleEnd)>0.998)
      )
    {
      return true;
    }
  }
  return false;
}

  double MFL::MinTime(std::vector<ana::hit> h){
    double ret = std::numeric_limits<double>::max();
    for(auto & _h : h)
    {
      ret=std::min(ret,_h.get_Time());
    }
    return ret;
  }

std::vector<ana::hit> MFL::GetHitsForTrack(art::Ptr<recob::Track> ptrack, art::Event const& e, art::InputTag trackModuleLabel)
{
  std::vector<ana::hit> ret;
  //moved to header.
  //protoana::ProtoDUNETrackUtils        trackUtil; //This should probably move somewhere else where it can initialize less often.
  const std::vector<const recob::Hit*> Hits = trackUtil.GetRecoTrackHits(*ptrack, e, trackModuleLabel.label());
  bool time_cut_pass = false;
  for(auto& hit : Hits)
  {
    double t = hit->PeakTime();
    if( t < mfl_conf.minhitpeakcut || t > mfl_conf.maxhitpeakcut )
    {
      continue;
      //break;
    }
    time_cut_pass=true;
    ret.emplace_back(hit->PeakTime(), hit->Channel(), hit->Integral(), hit->PeakAmplitude());
  }
  //Add HitTimesCuts here TODO
  if(time_cut_pass)
  {
    std::sort(ret.begin(), ret.end(), [](auto& a, auto& b)->bool{return a.get_Time()<b.get_Time();});
    //Get all hits for track using Protodune utils. Fill ret, sort ret, done.
    return ret;
  }else{
    ret.clear(); //return an empty vector.
    return ret;
  }
}


std::vector<ana::ohit> MFL::GetOpHitsForTrack(art::Ptr<recob::Track> ptrack, double T00, double ophitTrigT, std::vector<art::Ptr<recob::OpHit> > const& ophitList)
{
  std::vector<ana::ohit> ret;
  for(auto const& hit : ophitList)
  {
    double T0 = hit->PeakTime() - ophitTrigT;
    if( T0 < T00+mfl_conf.OpHitWindowMax && T0 > T00 - mfl_conf.OpHitWindowMin)
    {
      ret.emplace_back(hit->PeakTime(), hit->OpChannel(), hit->Width(), hit->Area(), hit->PE());
      std::cout<<"OPHIT PASSED\n";
    }else{
      std::cout<<"OPHIT FAILED\n";
    }
    std::sort(ret.begin(), ret.end(), [](auto const& a, auto const& b)->bool {return a.get_Time() < b.get_Time();});
  }
  return ret;
}

std::vector<ana::oflash> MFL::GetOpFlashesForTrack(art::Ptr<recob::Track> ptrack, double T00, double opflashTrigT, std::vector<art::Ptr<recob::OpFlash> > const& opflashList)
{
  std::vector<ana::oflash> ret;
  for(auto const& flash : opflashList)
  {
    double T0 = flash->Time() - opflashTrigT;
    if( T0 < T00+mfl_conf.OpFlashWindowMax && T0 > T00 - mfl_conf.OpFlashWindowMin)
    {
      double sum(0);
      auto tmp = flash->PEs();
      for(auto const& v : tmp)
      {
        sum+=v;
      }
      ret.emplace_back(flash->Time(), flash->TimeWidth(), sum, tmp, flash->FastToTotal());
    }
    std::sort(ret.begin(), ret.end(), [](auto const& a, auto const& b)->bool {return a.get_Time() < b.get_Time();});
  }
  return ret;
}


/*
   std::vector<ana::hit> MFL::CandidateCollectionHits(std::vector<ana::hit> & hits, art::FindManyP<recob::Hit, recob::TrackHitMeta> mp_thm, size_t trk_size)
   {
   std::vector<ana::hit> tmp;
   tmp.resize(hits.size());
   std::copy(hits.begin(), hits.end(), tmp.begin());
   for(size_t i=0; i<hits.size(); i++)
   {

   }
   for(size_t i=0; i<trk_size; i++)
   {

   }//WORKING HERE
   return hits;
   }
   */

//TODO: This function is incomplete, and should also be refactored. This is ugly, but I haven't finished making sense of Alena's methods here.
std::vector<std::vector<double>> MFL::HasMichel( std::vector<art::Ptr<recob::Track>> const& tracklist, art::FindManyP<recob::Hit, recob::TrackHitMeta> const& fmthm, TVector3 const& pos, TVector3 const& end, size_t const& trknum, double const& T00, art::FindManyP<recob::SpacePoint> const& fmsp, art::Handle<std::vector<recob::Track>> trackListHandle, std::vector<art::Ptr<recob::Hit>> const& hitlist, std::vector<art::Ptr<recob::Shower>> const& showerlist, art::FindManyP<recob::PFParticle> const& pfp_shwr_assn, art::FindManyP<anab::T0> const& shwr_t0_assn_v, detinfo::DetectorPropertiesData const& detprop )
{
  const int kMaxHits =  mfl_conf.kMaxHits;
  std::vector<std::vector<double>> ret;
  int kMaxCh = mfl_conf.kMaxCh;
  double _hitdist = mfl_conf.hitdist;
  double id(-99999);
  double dist(id), enddist(id), endtpcno(id);//, endhitkey(id),
  //       endwireno(id), endchno(id), endhitchrg(id), endpeaktime(id),
  //      endhitx(id), endhity(id), endhitz(id);
  std::vector<double> hits_key(kMaxCh,-999);
  std::vector<double> hits_charge(kMaxCh,-999);
  std::vector<double> hits_chno(kMaxCh,-999);
  std::vector<double> hits_wire(kMaxCh,-999);
  std::vector<double> hits_peakT(kMaxCh,-999);
  std::vector<double> hits_TPC(kMaxCh,-999);
  std::vector<double>  hits_xpos(kMaxCh,-999);
  std::vector<double>  hits_ypos(kMaxCh,-999);
  std::vector<double>  hits_zpos(kMaxCh,-999);
  int APAnegbound1=mfl_conf.APAnegbound1();
  int APAposbound1=mfl_conf.APAposbound1();
  int APAnegbound2=mfl_conf.APAnegbound2();
  int APAposbound2=mfl_conf.APAposbound2();
  //float trkbegx = pos.X();
  //float trkbegy = pos.Y();
  //float trkbegz = pos.Z();
  double trkstopx = end.X();
  double trkstopy = end.Y();
  double trkstopz = end.Z();
  double otherlongtrklen=mfl_conf.oltl;


  std::vector <double> trkHitsKey;
  for(size_t oo=0; oo<tracklist.size(); oo++)
  {
    art::Ptr<recob::Track> ptrack_l(trackListHandle, oo);
    const recob::Track& track_l = *ptrack_l;
    auto vhit=fmthm.at(oo);
    auto vmeta=fmthm.data(oo);
    if(track_l.Length()> otherlongtrklen )
    {
      for (size_t zz = 0; zz<vhit.size(); ++zz)
      {
        if(vhit[zz]->WireID().Plane == 2)
        {
          trkHitsKey.push_back(vhit[zz].key());
        }
      }//loop over hits
    }//if(track_l.Length()>75 )
  }//loop over tracks


  int trkcolhits{0};
  std::vector<double> longtrk_hitkey;

  if(fmthm.isValid())
  {
    auto vhit=fmthm.at(trknum);
    auto vmeta=fmthm.data(trknum);
    for(size_t i=0; i<vhit.size();++i)
    {
      //if wrong plane or if bad hit.
      if (vhit[i]->WireID().Plane!=2 ||
          vmeta[i]->Index() == std::numeric_limits<int>::max() ||
          vmeta[i]->Index()>=tracklist[i]->NumberTrajectoryPoints() ||
          !tracklist[i]->HasValidPoint(vmeta[i]->Index()) ){ continue;}

      TVector3 loc = tracklist[i]->LocationAtPoint<TVector3>(vmeta[i]->Index());
      longtrk_hitkey.push_back(vhit[i].key());
      if( !( loc.Z()!=0 && (loc.Z()<APAnegbound1 || loc.Z()>APAposbound1) && (loc.Z()<APAnegbound2 || loc.Z()>APAposbound2) ) ){continue;}
      hits_key.at(trkcolhits)       =   vhit[i].key();
      hits_charge.at(trkcolhits)    =   vhit[i]->Integral();
      hits_wire.at(trkcolhits)      =   vhit[i]->WireID().Wire;
      hits_peakT.at(trkcolhits)     =   vhit[i]->PeakTime();
      hits_TPC.at(trkcolhits)       =   vhit[i]->WireID().TPC;
      hits_chno.at(trkcolhits)      =   vhit[i]->Channel();
      hits_xpos.at(trkcolhits)      =   loc.X();
      hits_ypos.at(trkcolhits)      =   loc.Y();
      hits_zpos.at(trkcolhits)      =   loc.Z();
      dist=sqrt(pow(hits_xpos.at(trkcolhits)-trkstopx,2)+pow(hits_ypos.at(trkcolhits)-trkstopy,2)+pow(hits_zpos.at(trkcolhits)-trkstopz,2));
      trkcolhits++;
      if(abs(dist)>enddist){continue;}
      //enddist=dist;
      //endhitkey=vhit[i].key();
      //endwireno=vhit[i]->WireID().Wire;
      //endpeaktime=vhit[i]->PeakTime();
      //endtpcno=vhit[i]->WireID().TPC;
      //endchno=vhit[i]->Channel();
      //endhitchrg=vhit[i]->Integral();
      //endhitx=loc.X();
      //endhity=loc.Y();
      //endhitz=loc.Z();
    }

  }
  double dist0{999},dist1{999};
  if(trkcolhits==0){return ret;};//return empty vector. Don't know what to do if trkcolhits empty
  dist0=sqrt(pow(hits_xpos[0]-trkstopx,2)+pow(hits_ypos[0]-trkstopy,2)+pow(hits_zpos[0]-trkstopz,2));
  dist1=sqrt(pow(hits_xpos[trkcolhits-1]-trkstopx,2)+pow(hits_ypos[trkcolhits-1]-trkstopy,2)+pow(hits_zpos[trkcolhits-1]-trkstopz,2));
  int a(0),c(0),e(0);
  double b(0.0),d(0.0),f(0.0),g(0.0),h(0.0),r(0.0);
  if(abs(dist1)>abs(dist0))
  {
    for(int jj=0; jj<trkcolhits-1; jj++)
    {
      for (int kk = jj + 1; kk < trkcolhits; kk++)
      {
        a = hits_key[jj];
        hits_key[jj] = hits_key[kk];
        hits_key[kk] = a;

        b = hits_charge[jj];
        hits_charge[jj] = hits_charge[kk];
        hits_charge[kk] = b;

        c = hits_wire[jj];
        hits_wire[jj] = hits_wire[kk];
        hits_wire[kk] = c;

        d = hits_peakT[jj];
        hits_peakT[jj] = hits_peakT[kk];
        hits_peakT[kk] = d;

        e = hits_TPC[jj];
        hits_TPC[jj] = hits_TPC[kk];
        hits_TPC[kk] = e;

        f = hits_chno[jj];
        hits_chno[jj] = hits_chno[kk];
        hits_chno[kk] = f;

        g = hits_xpos[jj];
        hits_xpos[jj] = hits_xpos[kk];
        hits_xpos[kk] = g;

        h = hits_ypos[jj];
        hits_ypos[jj] = hits_ypos[kk];
        hits_ypos[kk] = h;

        r = hits_zpos[jj];
        hits_zpos[jj] = hits_zpos[kk];
        hits_zpos[kk] = r;
      }
    }
  }

  int nearhitct=0;
  std::vector<int> nearhits_key(kMaxHits,-999);
  std::vector<double> nearhits_peakT(kMaxHits, -999);
  std::vector<double> nearhits_charge (kMaxHits, -999);
  std::vector<double> nearhits_wire (kMaxHits, -999);
  std::vector<double> nearhits_chno (kMaxHits, -999);
  std::vector<double> nearhits_TPC  (kMaxHits, -999);
  std::vector<double> nearhits_plane(kMaxHits, -999);
  std::vector<double> nearhits_xpos (kMaxHits, -999);
  std::vector<double> nearhits_ypos (kMaxHits, -999);
  std::vector<double> nearhits_zpos (kMaxHits, -999);
  //  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  art::ServiceHandle<geo::Geometry const> geometry;

  for(size_t ll=0; ll<hitlist.size();++ll) //loop over all hits
  {
    //auto tracks = thass.at(hitlist[ll].key());


    if(hitlist[ll]->WireID().Plane==2)
    {
      //----Make sure that this hit does not belong to the candidate muon track----//
      int same_trk = 0;
      for(size_t mm=0;mm<longtrk_hitkey.size();mm++)
      {
        if(longtrk_hitkey.at(mm)==hitlist[ll].key()) { same_trk = 1; }
      }
      if ( same_trk==1 ) continue;

      //---Make sure that these hits don't belong to another long track-----//
      int long_trk = 0;
      for(size_t o1=0;o1<trkHitsKey.size();o1++)
      {
        if(trkHitsKey.at(o1)==hitlist[ll].key()) { long_trk = 1; }
      }
      if(long_trk==1) continue;

      double diffpeaktt0     = hitlist[ll]->PeakTime() - (T00/1000)*2; //calculate the peak time in ticks
      double allhitX = detprop.ConvertTicksToX(diffpeaktt0,hitlist[ll]->WireID().Plane,hitlist[ll]->WireID().TPC,hitlist[ll]->WireID().Cryostat);
      double Wirestart[3], Wireend[3];
      geometry->WireEndPoints(hitlist[ll]->WireID().Cryostat, hitlist[ll]->WireID().TPC, hitlist[ll]->WireID().Plane, hitlist[ll]->WireID().Wire, Wirestart,Wireend);
      double allhitZ = Wirestart[2];
      if(!fmsp.at(ll).size()) continue;
      double allhitY = fmsp.at(ll)[0]->XYZ()[1];
      double hitdist = sqrt(pow(abs(allhitX-trkstopx),2)+pow(abs(allhitZ-trkstopz),2));
      if( hitlist[ll]->WireID().TPC==endtpcno &&
          (allhitZ<APAnegbound1 || allhitZ>APAposbound1) &&
          (allhitZ<APAnegbound2 ||allhitZ>APAposbound2) &&
          hitdist <= _hitdist ) //  This is a hit that belongs to the michel cone.
      {
        nearhits_xpos.at(nearhitct)    = allhitX;
        nearhits_ypos.at(nearhitct)    = allhitY;
        nearhits_zpos.at(nearhitct)    = allhitZ;
        nearhitct++;
      }
    }//if(hitlist[ll]->WireID().Plane==2)
  }//for(size_t ll=0; ll<hitlist.size();++ll
  //double TPC_trigger_offset = 0.0;
  //auto const* detclock = lar::providerFrom<detinfo::DetectorClocksService>();
  //TPC_trigger_offset = detclock->TriggerOffsetTPC();

  double shwr_dist        = 99999, shwr_dist0 = 99999;
  //int    shwr_key         = -999;
  //int    shwr_ID          = -999;
  //double shwr_length      = -999;
  //double shwr_startx      = -999;
  //double shwr_starty      = -999;
  //double shwr_startz      = -999;
  //int    shwr_bestplane   = -999;
  //double shwr_startdcosx  = -999;
  //double shwr_startdcosy  = -999;
  //double shwr_startdcosz  = -999;
  //double shwr_openangle   = -999;


  for(size_t jjj=0; jjj<showerlist.size();++jjj)
  {
    //art::Ptr<recob::Shower> pshower(showerListHandle, jjj);
    const recob::Shower& shower = *(showerlist.at(jjj));

    std::vector<art::Ptr<recob::PFParticle>> pfpsh=pfp_shwr_assn.at(jjj);
    if(pfpsh.size())
    {

      //---Only consider T0 tagged shower-----//

      std::vector<art::Ptr<anab::T0>> t0sh=shwr_t0_assn_v.at(pfpsh[0].key());
      if(t0sh.size() )
      {
        TVector3 const& shwr_start = shower.ShowerStart();
        float shwr_startx1     = shwr_start.X();
        //TODO: WHy do we disregaurd Y? Can I turn it on in shwr_dist? JStock
        //NOTE*: We disregaurd Y in HitDist as well.
        //float shwr_starty1     = shwr_start.Y();
        float shwr_startz1     = shwr_start.Z();

        shwr_dist = sqrt(pow(shwr_startx1 - trkstopx,2) + pow(shwr_startz1 - trkstopz,2));

        if( shwr_dist<shwr_dist0) //This is the closest shower
        {
          shwr_dist0       = shwr_dist;
          std::cout<<"Closest Shower: "<<shwr_dist<<"\n";
          //  shwr_key         = showerlist[jjj].key();
          //  shwr_ID          = shower.ID();
          //  shwr_length      = shower.Length();
          //  shwr_startx      = shwr_startx1;
          //  shwr_starty      = shwr_starty1;
          //  shwr_startz      = shwr_startz1;
          //  shwr_bestplane   = shower.best_plane();
          //  TVector3 const& shwrdir_start = shower.Direction();
          //  shwr_startdcosx  = shwrdir_start.X();
          //  shwr_startdcosy  = shwrdir_start.Y();
          //  shwr_startdcosz  = shwrdir_start.Z();
          //  shwr_openangle   = shower.OpenAngle();
        }
      }
    }
  }

  return ret;
}

//TODO GetOpFlash
/*bool MFL::HitTimesCut(art::Ptr<recob::Track> const& ptrack, art::Event const& e, std::string label)
  {
  const std::vector<const recob::Hit*> Hits =
  trackUtil.GetRecoTrackHits(*ptrack, e, label);
  }*/


