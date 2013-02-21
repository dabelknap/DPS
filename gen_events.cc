/**
 * @author D. Austin Belknap
 * @file gen_events.cc
 *
 * The program simulates the HZZ4L Z-candidate
 * selection algorithm and stores the event
 * information in a ROOT TTree
 */

#include <iostream>
#include <cmath>
#include <algorithm>

#include "Pythia.h"
#include "TTree.h"
#include "TFile.h"

using namespace Pythia8;
using std::cout;
using std::endl;

const double ZMASS = 91.188;

/**
 * Allows the sort and next_permutation functions
 * to compare leptons by pT
 */
bool comp_leptons( Particle a, Particle b )
{
    return ( a.pT() > b.pT() );
}


int main( int argc, char* argv[] )
{
    if (argc < 3)
    {
        cout << "Usage: gen_events [config_file.cmnd] [out_file.root]" << endl;
        exit(1);
    }

    Pythia pythia;
    pythia.readFile( argv[1] );

    int nEvents = pythia.mode("Main:numberOfEvents");

    double z1l1pt, z1l1eta, z1l1phi;
    double z1l2pt, z1l2eta, z1l2phi;
    double z2l1pt, z2l1eta, z2l1phi;
    double z2l2pt, z2l2eta, z2l2phi;

    int z1l1ID, z1l2ID, z2l1ID, z2l2ID;

    double z1mass, z1pt, z1eta, z1phi;
    double z2mass, z2pt, z2eta, z2phi;

    double mass, genMass;

    TFile *f = new TFile( argv[2], "RECREATE");
    TTree *t = new TTree("selectedEvents","selectedEvents");


    // Book TTree Branches

    t->Branch("nEvents", &nEvents, "nEvents/I");

    t->Branch("mass", &mass, "mass/D");
    t->Branch("genMass", &mass, "genMass/D");
    t->Branch("z1mass", &z1mass, "z1mass/D");
    t->Branch("z2mass", &z2mass, "z2mass/D");

    t->Branch("z1pt", &z1pt, "z1pt/D");
    t->Branch("z1eta", &z1eta, "z1eta/D");
    t->Branch("z1phi", &z1phi, "z1phi/D");

    t->Branch("z2pt", &z2pt, "z2pt/D");
    t->Branch("z2eta", &z2eta, "z2eta/D");
    t->Branch("z2phi", &z2phi, "z2phi/D");

    t->Branch("z1l1pt", &z1l1pt, "z1l1pt/D");
    t->Branch("z1l1eta", &z1l1eta, "z1l1eta/D");
    t->Branch("z1l1phi", &z1l1phi, "z1l1phi/D");
    t->Branch("z1l1ID", &z1l1ID, "z1l1ID/I");

    t->Branch("z1l2pt", &z1l2pt, "z1l2pt/D");
    t->Branch("z1l2eta", &z1l2eta, "z1l2eta/D");
    t->Branch("z1l2phi", &z1l2phi, "z1l2phi/D");
    t->Branch("z1l2ID", &z1l2ID, "z1l2ID/I");

    t->Branch("z2l1pt", &z2l1pt, "z2l1pt/D");
    t->Branch("z2l1eta", &z2l1eta, "z2l1eta/D");
    t->Branch("z2l1phi", &z2l1phi, "z2l1phi/D");
    t->Branch("z2l1ID", &z2l1ID, "z2l1ID/I");

    t->Branch("z2l2pt", &z2l2pt, "z2l2pt/D");
    t->Branch("z2l2eta", &z2l2eta, "z2l2eta/D");
    t->Branch("z2l2phi", &z2l2phi, "z2l2phi/D");
    t->Branch("z2l2ID", &z2l2ID, "z2l2ID/I");

    pythia.init();

    int passed = 0;

    for ( int i = 0; i < nEvents; i++ )
    {
        pythia.next();

        std::vector<Particle> leptons;
        std::vector<Particle> final_leptons (4);

        genMass = 0;

        for ( int j = 0; j < pythia.event.size(); ++j )
        {
            Particle particle = pythia.event[j];

            if (particle.id() == 25)
                genMass = particle.mass();

            bool muon_pass     = abs(particle.id()) == 13 && particle.pT() > 5 && fabs(particle.eta()) < 2.4 && particle.isFinal();
            bool electron_pass = abs(particle.id()) == 11 && particle.pT() > 7 && fabs(particle.eta()) < 2.5 && particle.isFinal();

            if ( muon_pass || electron_pass )
                leptons.push_back( particle );
        }

        if ( leptons.size() < 4 )
            continue;

        std::sort( leptons.begin(), leptons.end(), comp_leptons );

        if ( !(leptons.at(0).pT() > 20 && leptons.at(1).pT() > 10 ) )
            continue;

        double prev_z1mass = 0;
        double prev_z2l1pt = 0;
        double prev_z2l2pt = 0;

        // permute over different arrangements of the final-state leptons in the event
        do
        {
            // match OS SF leptons
            bool pass_OSSF    = leptons.at(0).id() == -leptons.at(1).id() && leptons.at(2).id() == -leptons.at(3).id();

            // ensure leptons are in descending order per Z candidate
            bool pass_pTorder = leptons.at(0).pT() > leptons.at(1).pT() && leptons.at(2).pT() > leptons.at(3).pT();

            if ( ! (pass_OSSF && pass_pTorder) )
                continue;

            // Z1 is closes to nominal Z mass
            z1mass = m(leptons.at(0).p(),leptons.at(1).p());
            if ( fabs(z1mass-ZMASS) < fabs(prev_z1mass-ZMASS) )
            {
                prev_z1mass = z1mass;
                std::copy( leptons.begin(), leptons.begin() + 4, final_leptons.begin() );
            }

            // Z2 is made from highest pt leptons
            z2l1pt = leptons.at(2).pT();
            z2l2pt = leptons.at(3).pT();
            if ( z1mass == prev_z1mass && z2l1pt > prev_z2l1pt && z2l2pt > prev_z2l2pt )
            {
                prev_z2l1pt = z2l1pt;
                prev_z2l2pt = z2l2pt;
                std::copy( leptons.begin(), leptons.begin() + 4, final_leptons.begin() );
            }
        }
        while ( std::next_permutation( leptons.begin(), leptons.end(), comp_leptons ) );

        // continue only if four final leptons have been selected
        if ( final_leptons.size() < 4 )
            continue;

        z1l1pt = final_leptons.at(0).pT();
        z1l2pt = final_leptons.at(1).pT();
        z2l1pt = final_leptons.at(2).pT();
        z2l2pt = final_leptons.at(3).pT();

        z1l1eta = final_leptons.at(0).eta();
        z1l2eta = final_leptons.at(1).eta();
        z2l1eta = final_leptons.at(2).eta();
        z2l2eta = final_leptons.at(3).eta();

        z1l1phi = final_leptons.at(0).phi();
        z1l2phi = final_leptons.at(1).phi();
        z2l1phi = final_leptons.at(2).phi();
        z2l2phi = final_leptons.at(3).phi();

        z1l1ID = final_leptons.at(0).id();
        z1l2ID = final_leptons.at(1).id();
        z2l1ID = final_leptons.at(2).id();
        z2l2ID = final_leptons.at(3).id();

        Vec4 z1p4 = final_leptons.at(0).p() + final_leptons.at(1).p();
        Vec4 z2p4 = final_leptons.at(2).p() + final_leptons.at(3).p();

        z1mass = z1p4.mCalc();
        z1pt   = z1p4.pT();
        z1eta  = -log( tan(z1p4.theta()/2.0) );
        z1phi  = z1p4.phi();

        z2mass = z2p4.mCalc();
        z2pt   = z2p4.pT();
        z2eta  = -log( tan(z2p4.theta()/2.0) );
        z2phi  = z2p4.phi();

        mass   = (z1p4 + z2p4).mCalc();

        t->Fill();

    }
    t->Write();
    f->Close();

    pythia.stat();

    return 0;
}
