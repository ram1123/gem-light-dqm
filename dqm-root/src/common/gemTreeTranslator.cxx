#ifndef DEBUG
#define DEBUG 1
#endif
#include <TChain.h>
#include <TProofOutputFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <sstream>
#include <vector>
#include <TSystem.h>
#include <bitset>

#include "AMC13_histogram.cxx"



unsigned int countBits(uint64_t n)
{
    unsigned int count = 0;
        while (n)
        {
            count += n & 1;
            n >>= 1;
        }

return count;
}
class gemTreeTranslator: public TSelector {
    public :
        TTree          *fChain;   //!pointer to the analyzed TTree or TChain
        TTreeReader     fReader;  //!the tree reader

        // Declaration of leaf types
        Event           *GEMEvents;
        // List of branches
        TBranch        *b_GEMEvents;   //!
        gemTreeTranslator(TTree * /*tree*/ =0) : fChain(0) { }
        //gemTreeTranslator(TTree * /*tree*/ =0) : { }
        virtual ~gemTreeTranslator() { }
        virtual Int_t   Version() const { return 2; }
        virtual void    Begin(TTree *tree);
        virtual void    SlaveBegin(TTree *tree);
        virtual void    Init(TTree *tree);
        virtual Bool_t  Notify();
        virtual Bool_t  Process(Long64_t entry);
        virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
        virtual void    SetOption(const char *option) { fOption = option; }
        virtual void    SetObject(TObject *obj) { fObject = obj; }
        virtual void    SetInputList(TList *input) { fInput = input; }
        virtual TList  *GetOutputList() const { return fOutput; }
        virtual void    SlaveTerminate();
        virtual void    Terminate();

        vector<AMC13Event> v_amc13;    ///<Vector of AMC13Event
        vector<AMCdata> v_amc;         ///<Vector of AMCdata
        vector<GEBdata> v_geb;         ///<Vector of GEBdata
        vector<VFATdata> v_vfat;       ///Vector of VFATdata

        //AMC_histogram * v_amcH;        ///<Vector of AMC_histogram
        //GEB_histogram * v_gebH;        ///<Vector of GEB_histogram
        //VFAT_histogram * v_vfatH;      ///<Vector of VFAT_histogram

        //BRS
        TTree* treeTranslator; //Should make this private
        std::vector<TTree*> tree_list;
        //int amcID;
        //int amcBID;
        //int gebID;
        //int infoAMCVEC;
        //int infoGEBVEC;
        //int infoVFATVEC;

        int calPhase;
        //int Dly;
        int l1aTime;
        int latency;
        int link;
        //int pDel;
        //int mode;
        //int mspl;
        //int Nev;
        int Nhits;
        //int trimDAC;
        //int utime;
        //int vcal;
        //int vfatCH;
        //int vfatID;
        int vfatN;
        int vth;
        int vth1;
        //int vth2;
        //float ztrim;
        //int shelf;
        int slot;



        //int VFATMap[12][12][24];

        //  AMC13_histogram * m_amc13H;
        //  AMC_histogram * m_amcH;
        //  GEB_histogram * m_gebH;
        //  VFAT_histogram * m_vfatH;

        //int m_deltaV;
        //int m_Latency;
        //long long int m_OrbitNumber;
        //long long int m_RelOrbitNumber;
        //TDirectory* onlineHistsDir;

        TFile *fFile;
        TProofOutputFile *fProofFile;

        ClassDef(gemTreeTranslator,2);
};

void gemTreeTranslator::Begin(TTree * /*tree*/)
{
    // The Begin() function is called at the start of the query.
    // When running with PROOF Begin() is only called on the client.
    // The tree argument is deprecated (on PROOF 0 is passed).

    if (DEBUG) std::cout << "MASTER BEGIN"<< std::endl;
    TString option = GetOption();
}

void gemTreeTranslator::SlaveBegin(TTree * /*tree*/)
{
    // The SlaveBegin() function is called after the Begin() function.
    // When running with PROOF SlaveBegin() is called on each slave server.
    // The tree argument is deprecated (on PROOF 0 is passed).

    if (DEBUG) std::cout << "SLAVE BEGIN"<< std::endl;
    TString option = GetOption();
    gSystem->Load("libEvent.so");

    // The file for merging
    fProofFile = new TProofOutputFile("SimpleFile.root", TProofOutputFile::kMerge);
    TNamed *out = (TNamed *) fInput->FindObject("PROOF_OUTPUTFILE"); //analyzed.root
    if (out) fProofFile->SetOutputFileName(out->GetTitle());
    TDirectory *savedir = gDirectory;
    fFile = fProofFile->OpenFile("RECREATE");
    if (fFile && fFile->IsZombie()) SafeDelete(fFile);
    savedir->cd();

    // Cannot continue
    if (!fFile) {
        TString amsg = TString::Format("ProofSimpleFile::SlaveBegin: could not create '%s':"
                " instance is invalid!", fProofFile->GetName());
        Abort(amsg, kAbortProcess);
        return;
    }
    fFile->cd();

    // BRS: Lets make the Tfile here:
    int numOfAMC=3;
    int numOfGEB=5;
    int numOfVFAT=24;

    tree_list.clear(); // clear the vector having list of trees
    for (auto amc_count=0; amc_count<numOfAMC; ++amc_count){
        for (auto geb_count=0; geb_count<numOfGEB; ++geb_count){
            //for (auto vfat_count=0; vfat_count<numOfVFAT; ++vfat_count){
                TString treeName="gemTree_"+std::to_string(amc_count)+"_"+std::to_string(geb_count);
                new TProofOutputFile(treeName, "M");
                treeTranslator    = new TTree(treeName,"Tree holding gem info");
                auto branchcalPhase  = treeTranslator->Branch("calPhase",&calPhase,"calPhase/I"); //amcID
                //auto branchDly  = treeTranslator->Branch("Dly",&Dly,"Dly/I"); //amcID
                auto branchl1aTime  = treeTranslator->Branch("l1aTime",&l1aTime,"l1aTime/I"); //amcID
                auto branchlatency  = treeTranslator->Branch("latency",&latency,"latency/I"); //amcID
                auto branchlink  = treeTranslator->Branch("link",&link,"link/I"); //amcID
                //auto branchpDel  = treeTranslator->Branch("pDel",&pDel,"pDel/I"); //amcID
                //auto branchmode  = treeTranslator->Branch("mode",&mode,"mode/I"); //amcID
                //auto branchmspl  = treeTranslator->Branch("mspl",&mspl,"mspl/I"); //amcID
                //auto branchNev  = treeTranslator->Branch("Nev",&Nev,"Nev/I"); //amcID
                auto branchNhits  = treeTranslator->Branch("Nhits",&Nhits,"Nhits/I"); //amcID
                //auto branchtrimDAC  = treeTranslator->Branch("trimDAC",&trimDAC,"trimDAC/I"); //amcID
                //auto branchutime  = treeTranslator->Branch("utime",&utime,"utime/I"); //amcID
                //auto branchvcal  = treeTranslator->Branch("vcal",&vcal,"vcal/I"); //amcID
                //auto branchvfatCH  = treeTranslator->Branch("vfatCH",&vfatCH,"vfatCH/I"); //amcID
                //auto branchvfatID  = treeTranslator->Branch("vfatID",&vfatID,"vfatID/I"); //amcID
                auto branchvfatN  = treeTranslator->Branch("vfatN",&vfatN,"vfatN/I"); //amcID
                auto branchvth  = treeTranslator->Branch("vth",&vth,"vth/I"); //amcID
                auto branchvth1  = treeTranslator->Branch("vth1",&vth1,"vth1/I"); //amcID
                //auto branchvth2  = treeTranslator->Branch("vth2",&vth2,"vth2/I"); //amcID
                //auto branchztrim  = treeTranslator->Branch("ztrim",&ztrim,"ztrim/F"); //amcID


                //auto branchshelf  = treeTranslator->Branch("shelf",&shelf,"shelf/I"); //amcID
                auto branchslot  = treeTranslator->Branch("slot",&slot,"slot/I"); //amcID

                //auto branchAMCID  = treeTranslator->Branch("amcID",&amcID,"AMC ID/I"); //amcID
                //auto branchGEBID  = treeTranslator->Branch("gebID",&gebID,"GEB ID/I");
                tree_list.push_back(treeTranslator);
            //}  // End of VFAT loop
        } // End of GEB loop
    }  // End of AMC loop

    //fFile->Add(treeTranslator);
    for (auto obj : tree_list) {
        fFile->Add(obj);
    }

    if (fChain){
        fChain->GetEntry(0);
        v_amc13 = GEMEvents->amc13s();
        if (DEBUG) std::cout << "SlaveBegin: v_amc13 = "<< v_amc13.size() << std::endl;
    }
    else{
        Error("Begin", "no fChain!");
    }
    ///// AMCdata::BID() or AMCdata::AMCNum() //AMC ID..
    ///// VFATdata::ChipID() vfatID
    ///// shelf and slot?

    if (DEBUG) std::cout << "SLAVE END"<< std::endl;
}

//!Fills the histograms that were book from bookAllHistograms
Bool_t gemTreeTranslator::Process(Long64_t entry)
{
    //fReader.SetEntry(entry);
    fChain->GetEntry(entry);
    if (DEBUG) std::cout<<"########################################"<<std::endl;
    if (DEBUG) std::cout<<"\n\tEvent loop\n"<<std::endl;
    if (DEBUG) std::cout<<"########################################"<<std::endl;
    int amc13_count=0;    //counter through AMC13s
    int amc_count=0;      //counter through AMCs
    int geb_count=0;      //counter through GEBs
    int vfat_count=0;      //counter through VFATs

    /* local variables that will be filled to the tree  */
    int local_amcID = -1;
    int local_calPhase = -1;
    int local_Dly = -1;
    int local_l1aTime = -1;
    int local_latency = -1;
    int local_link = -1;
    int local_shelf = -1;
    int local_slot = -1;
    //int local_vfatID = -1;
    int local_vfatN = -1;
    int local_vth = -1;
    int local_vth1 = -1;
    //int local_gebID = -1;
    int local_Nev = 0;
    int local_Nhits = 0;
    /* END local variable*/

    v_amc13 = GEMEvents->amc13s();
    if (DEBUG) cout << "Get a vector of AMC13 having size " << v_amc13.size() << endl;
    /* LOOP THROUGH AMC13s */
    for(auto amc13 = v_amc13.begin(); amc13!=v_amc13.end(); amc13++) {
        v_amc = amc13->amcs();
        if (DEBUG) cout << "Found AMC having size " << v_amc.size() << endl;
        //vector<unsigned short> test = amc13->BoardID();
        //if (DEBUG) cout << "BoardID = "<< (amc13->BoardID()).size() << endl;
        //local_shelf = amc13->BoardID();
        /* LOOP THROUGH AMCs */
        for(auto amc=v_amc.begin(); amc!=v_amc.end(); amc++){
            v_geb = amc->gebs();
            if (DEBUG) cout << "Found GEB having size "<< v_geb.size() << endl;
            local_slot =amc->AMCnum();
            local_l1aTime = amc->L1AT();
            int m_RunType = amc->Rtype();
            if (m_RunType) {
                local_latency = amc->Param1(); //BRS: NO IDEA
                local_vth = amc->Param2();
                local_vth1 = amc->Param3();
            }

            geb_count=0;
            //if (DEBUG) std::cout<<"amc:"<< amc_count << ": local_slot = "<<local_slot<<std::endl;
            //if (DEBUG) std::cout<<"amc:"<< amc_count << ": local_l1aTime = "<<local_l1aTime<<std::endl;
            //if (DEBUG) std::cout<<"amc:"<< amc_count << ": local_m_RunType= "<<m_RunType<<std::endl;
            //if (DEBUG) std::cout<<"amc:"<< amc_count << ": local_latency = "<<local_latency<<std::endl;
            //if (DEBUG) std::cout<<"amc:"<< amc_count << ": local_vth = "<<local_vth<<std::endl;
            //if (DEBUG) std::cout<<"amc:"<< amc_count << ": local_vth1 = "<<local_vth1<<std::endl;
            //if (DEBUG) std::cout<<"amc:  = "<<<<std::endl;
            /* LOOP THROUGH GEBs */
            for(auto g=v_geb.begin(); g!=v_geb.end();g++) {
                v_vfat = g->vfats();
                local_link = g->InputID();
                /* LOOP THROUGH VFATs */
                vfat_count = 0;
                if (DEBUG) std::cout<<"-----"<<endl;
                if (DEBUG) std::cout<<"GEB:" << geb_count << ": local_link = "<<local_link<<std::endl;
                if (DEBUG) std::cout << "Found VFATs having size "<< v_vfat.size() << std::endl;
                for(auto v=v_vfat.begin(); v!=v_vfat.end();v++) {
                    //if (DEBUG) std::cout << "Found VFATs " << std::endl;
                    TString treeName="gemTree_"+std::to_string(amc_count)+"_"+std::to_string(geb_count);
                    local_amcID = amc->AMCnum();
                    local_calPhase = v->isBlockGood();
                    //local_vfatID = v->chipID();
                    local_vfatN = v->Pos();
                    uint64_t binChannel_0_63 = v->lsData();
                    uint64_t binChannel_64_127 = v->msData();
                    local_Nev  = v->EC();
                    local_Nhits += countBits(binChannel_0_63) + countBits(binChannel_64_127);
                    
                    //if (DEBUG) std::cout<< "\n\n==========================\n\n" << std::endl;
                    //if (DEBUG) std::cout<<"VFAT:" << vfat_count << ": treeName = "<<treeName<<std::endl;
                    //if (DEBUG) std::cout<<"VFAT:" << vfat_count << ": local_amcID = "<<local_amcID<<std::endl;
                    //if (DEBUG) std::cout<<"VFAT:" << vfat_count << ": local_calPhase = "<<local_calPhase<<std::endl;
                    //if (DEBUG) std::cout<<"VFAT:" << vfat_count << ": local_Pos = "<<local_vfatN<<std::endl;
                    //if (DEBUG) std::cout<<"VFAT:" << vfat_count << ": local_SlotNumber = "<<v->SlotNumber()<<std::endl;
                    //if (DEBUG) std::cout<<"VFAT:" << vfat_count << ": countBits(binChannel_0_63) = "<<countBits(binChannel_0_63)<<std::endl;
                    //if (DEBUG) std::cout<<"VFAT:" << vfat_count << ": countBits(binChannel_64_127) = "<<countBits(binChannel_64_127)<<std::endl;
                    //if (binChannel_0_63 != 0 || binChannel_64_127 != 0) {
                    if (DEBUG) std::cout<<"VFAT:" << vfat_count << ": local_Pos = "<<local_vfatN<<std::endl;
                    //if (DEBUG) std::cout<<"VFAT:" << vfat_count <<  ": lsData = " << binChannel_0_63 << std::endl;
                    //if (DEBUG) std::cout<<"VFAT:" << vfat_count <<  ": msData = " << binChannel_64_127 << std::endl;
                    //}
                    //if (DEBUG) std::cout<< "\n\n==========================\n\n" << std::endl;
                    //m_Latency = amc->Param1(); //BRS: NO IDEA

                    for (auto obj : tree_list) {
                        if (obj->GetName() == treeName) {
                            static_cast<TTree*>(obj)->Fill();
                            //if (DEBUG) cout << "Object name : " << obj->GetName() << endl;
                        }
                    }
                    //
                    vfat_count++;
                } /* END VFAT LOOP */
                geb_count++;
            } /* END GEB LOOP */
            amc_count++;
        } /* END AMC LOOP */
        amc13_count++;
    } /* END AMC13 LOOP */
    /*  Save Branches */
    calPhase = local_calPhase;
    l1aTime = local_l1aTime;
    latency = local_latency;
    link = local_link;
    //Nev = local_Nev;
    Nhits = local_Nhits;
    //vfatID = local_vfatID;
    vfatN = local_vfatN;
    vth = local_vth;
    vth1 = local_vth1;
    //shelf = local_shelf;
    slot = local_slot;

    //gebID = local_gebID;
    //amcID = local_amcID;

    treeTranslator->Fill();
    return kTRUE;
}

void gemTreeTranslator::SlaveTerminate()
{
    // The SlaveTerminate() function is called after all entries or objects
    // have been processed. When running with PROOF SlaveTerminate() is called
    // on each slave server.
    //
    // Should we fill summary canvases here?
    //
    //for (unsigned int i = 0; i < 12; i++){
    //    if (auto a = m_amc13H->amcsH(i)){
    //        for (unsigned int j = 0; j < 12; j++){
    //            if (auto g = a->gebsH(j)){
    //                g->fillSummaryCanvases();
    //                for (unsigned int k = 0; k < 24; k++){
    //                    if (auto v = g->vfatsH(k)){
    //                        v->fillWarnings();
    //                    }
    //                }
    //            }
    //        }
    //    }
    //}
    TDirectory *savedir = gDirectory;
    fFile->cd();
    fFile->Write();
    fOutput->Add(fProofFile);
    fFile->Close();
}
void gemTreeTranslator::Terminate()
{
    // The Terminate() function is the last function to be called during
    // a query. It always runs on the client, it can be used to present
    // the results graphically or save the results to file.
    Printf("Testing if gemTree was created successfully ....");
    TTree* gemTree = dynamic_cast<TTree*>(fOutput->FindObject("gemTree"));
    if (gemTree){
        Printf("gemTree exists!");
    }

    //  TTree* gemTreeB = dynamic_cast<TTree*>(fProofFile->FindObject("gemTree"));
    //if (treeTranslator){
    //Printf("gemTree exists! in memory at Terminate");
    //   treeTranslator->Print();
    //   treeTranslator->Show(0);
    //}
    //for (auto obj : tree_list)
    for(unsigned int i=0; i<tree_list.size(); i++) {
        treeTranslator = (TTree*)tree_list.at(i);
        treeTranslator->Write("", TObject::kOverwrite);
        delete treeTranslator;
        treeTranslator = 0;
    }

    if ((fProofFile = dynamic_cast<TProofOutputFile*>(fOutput->FindObject("SimpleFile.root")))) {
        TString outputFile(fProofFile->GetOutputFileName());
        TString outputName(fProofFile->GetName());
        Printf("outputFile: %s", outputFile.Data());
    }
}
//int gemTreeTranslator::slotFromMap(int a, int g, int cid)
//{
//  int res = -1;
//  for (int i=0; i<24; i++){
//    if (VFATMap[a][g][i] == cid) {res = i;}
//  }
//  return res;
//}
void gemTreeTranslator::Init(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;
    fChain->SetMakeClass(0);

    GEMEvents = new Event();

    fChain->SetBranchAddress("GEMEvents", &GEMEvents, &b_GEMEvents);
    v_amc13 = GEMEvents->amc13s();
    if (DEBUG) std::cout << "v_amc13 size = " << v_amc13.size() << std::endl;
}

Bool_t gemTreeTranslator::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}
