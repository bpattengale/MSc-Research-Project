// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/Beam.hh"
/// #include "Rivet/Projections/DressedLeptons.hh"
/// #include "Rivet/Projections/FinalState.hh"
/// #include "Rivet/Projections/MissingMomentum.hh"
/// #include "Rivet/Projections/DirectFinalState.hh"

namespace Rivet
{

      /// Analysis for 2007 Dijet Photoproduction Paper of HERA e-P and e+P collisions from 1998-2000
      class ZEUS_2007_I753991 : public Analysis
      {
      public:
            /// Constructor
            RIVET_DEFAULT_ANALYSIS_CTOR(ZEUS_2007_I753991);

            /// Book histograms and initialise projections before the run
            void init()
            {

                  // Initialise and register projections
                  FinalState fs;
                  declare(FastJets(fs, fastjet::JetAlgorithm::kt_algorithm, fastjet::RecombinationScheme::Et_scheme, 1.0), "Jets");
                  declare(DISKinematics(), "Kinematics");

                  // Histograms
                  // Table 2
                  book(_h_etbar[1], 1, 1, 1);

                  // Table 3
                  book(_h_etbar[0], 2, 1, 1);

                  // Table 4
                  book(_h_etjet1[1], 3, 1, 1);

                  // Table 5
                  book(_h_etjet1[0], 4, 1, 1);

                  // Table 6
                  book(_h_etabar[1], 5, 1, 1);

                  // Table 7
                  book(_h_etabar[0], 6, 1, 1);

                  // Table 8
                  book(_h_xpobs[1], 7, 1, 1);

                  // Table 9
                  book(_h_xpobs[0], 8, 1, 1);

                  // Table 10
                  book(_h_deltaphi[1], 9, 1, 1);

                  // Table 11
                  book(_h_deltaphi[0], 10, 1, 1);

                  // Table 12
                  book(_h_xpobs2, 11, 1, 1);

                  // Table 13
                  book(_h_xpobs3, 12, 1, 1);

                  // Table 14
                  book(_h_xpobs4, 13, 1, 1);

                  // Table 15
                  book(_h_xpobs5, 14, 1, 1);

                  // Table 16
                  book(_h_xpobs6, 15, 1, 1);

                  // Table 17
                  book(_h_xpobs7, 16, 1, 1);

                  // Table 18
                  book(_h_xpobs8, 17, 1, 1);

                  // Table 19
                  book(_h_xpobs9, 18, 1, 1);

                  // Table 20
                  book(_h_xyobs, 19, 1, 1);
            }

            /// Perform the per-event analysis
            void analyze(const Event &event)
            {

                  // Determine kinematics
                  const DISKinematics &kin = apply<DISKinematics>(event, "Kinematics");
                  if (kin.failed())
                        vetoEvent;
                  const int orientation = kin.orientation();

                  // Q2 cut and inelasticity cut
                  if (kin.Q2() > 1 * GeV2)
                        vetoEvent;
                  if (!inRange(kin.y(), 0.2, 0.85))
                        vetoEvent;

                  // Jet Selection
                  const Jets jets = apply<FastJets>(event, "Jets")
                                        .jets(Cuts::Et > 15 * GeV && Cuts::etaIn(-1 * orientation, 3 * orientation), cmpMomByEt);
                  MSG_DEBUG("Jet Multiplicity = " << jets.size());
                  if (jets.size() < 2) vetoEvent;
                  const Jet &j1 = jets[0];
                  const Jet &j2 = jets[1];
                  if (j1.Et() < 20 * GeV) vetoEvent;
                  const double eta1 = orientation * j1.eta();
                  const double eta2 = orientation * j2.eta();
                  if ((eta1 > 2.5) && (eta2 > 2.5)) vetoEvent;

                  // Jet eta bar and E bar calculation
                  const double etabar = (eta1 + eta2) / 2;
                  const double ebar = (j1.Et() + j2.Et()) / 2;

                  // Calculate x_y^obs and x_p^obs
                  const double xyobs = (j1.Et() * exp(-eta1) + j2.Et() * exp(-eta2)) / (2 * kin.y() * kin.beamLepton().E());
                  const size_t i_xyobs = (xyobs < 0.75) ? 0 : 1;
                  const double xpobs = (j1.Et() * exp(eta1) + j2.Et() * exp(eta2)) / (2 * kin.beamHadron().E());

                  // Delta Phi Calculation
                  const double phi1 = j1.phi();
                  const double phi2 = j2.phi();
                  double dPhi = phi1 - phi2;
                  if (dPhi > 3.14159265359)
                  {
                        dPhi = dPhi - 2*3.14159265359;
                  }
                  else if (dPhi < -3.14159265359)
                  {
                        dPhi = dPhi + 2*3.14159265359;
                  }
                  const double deltaPhi = abs(dPhi);

                  // Fill Histograms
                  _h_etbar[i_xyobs]->fill(ebar / GeV);
                  _h_etjet1[i_xyobs]->fill(j1.Et() / GeV);
                  _h_etabar[i_xyobs]->fill(etabar);
                  _h_xpobs[i_xyobs]->fill(xpobs);
                  _h_deltaphi[i_xyobs]->fill(deltaPhi);

                  // Symmetrize histograms with different eta regions
                  for(size_t isel = 0; isel < 2; ++isel)
                  {
                        double etaJet1 = (isel == 0) ? orientation*j1.eta() : orientation*j2.eta();
                        double etaJet2 = (isel == 0) ? orientation*j2.eta() : orientation*j1.eta();

                        if (i_xyobs > 0.75)
                        {
                              if (inRange(etaJet1, 0, 1) && inRange(etaJet2, 2, 3))
                              {
                                    if (j1.Et() > 25*GeV)
                                    {
                                          _h_xpobs2->fill(xpobs);
                                    }
                                    if (j1.Et() > 20*GeV)
                                    {
                                          _h_xpobs3->fill(xpobs);
                                    }
                              }
                              else if (inRange(etaJet1, -1, 0) && inRange(etaJet2, 0, 1))
                              {
                                    if (j1.Et() > 20*GeV)
                                    {
                                          _h_xpobs5->fill(xpobs);
                                    }
                              }
                        }
                        else
                        {
                              if (inRange(etaJet1, 2, 2.5) && inRange(etaJet2, 2, 3))
                              {
                                    if (j1.Et() > 20*GeV)
                                    {
                                          _h_xpobs6->fill(xpobs);
                                    }
                              }
                              else if (inRange(etaJet1, 1, 2) && inRange(etaJet2, 2, 3))
                              {
                                    if (j1.Et() > 20*GeV)
                                    {
                                          _h_xpobs8->fill(xpobs);
                                    }
                                    if (j1.Et() > 25*GeV)
                                    {
                                          _h_xpobs9->fill(xpobs);
                                    }
                              }
                        }
            	}

                  if (i_xyobs > 0.75)
                  {
                        if(inRange(eta1, 1, 2) && inRange(eta2, 1, 2))
                        {
                              if(j1.Et() > 30*GeV)
                              {
                                    _h_xpobs4->fill(xpobs);
                              }
                        }
                  }
                  else
                  {
                        if(inRange(eta1, 1, 2) && inRange(eta2, 1, 2))
                        {
                              if(j1.Et() > 25*GeV)
                              {
                                    _h_xpobs7->fill(xpobs);
                              }
                        }
                  }

                  _h_xyobs->fill(xyobs);
            }

            /// Normalise histograms etc., after the run
            void finalize()
            {
                  const double sf = crossSection() / picobarn / sumOfWeights();
                  for (auto &h : _h_etbar)
                        scale(h, sf);
                  for (auto &h : _h_etjet1)
                        scale(h, sf);
                  for (auto &h : _h_etabar)
                        scale(h, sf);
                  for (auto &h : _h_xpobs)
                        scale(h, sf);
                  for (auto &h : _h_deltaphi)
                        scale(h, sf);

                  scale(_h_xpobs2, sf);
                  scale(_h_xpobs3, sf);
                  scale(_h_xpobs4, sf);
                  scale(_h_xpobs5, sf);
                  scale(_h_xpobs6, sf);
                  scale(_h_xpobs7, sf);
                  scale(_h_xpobs8, sf);
                  scale(_h_xpobs9, sf);
                  scale(_h_xyobs, sf);
            }

            /// @}

      private:
            /// @name Histograms
            /// @{
            Histo1DPtr _h_etbar[2], _h_etjet1[2], _h_etabar[2], _h_xpobs[2], _h_deltaphi[2], _h_xpobs2, _h_xpobs3, _h_xpobs4, _h_xpobs5, _h_xpobs6, _h_xpobs7, _h_xpobs8, _h_xpobs9, _h_xyobs;
            /// @}
      };

      RIVET_DECLARE_PLUGIN(ZEUS_2007_I753991);

}
