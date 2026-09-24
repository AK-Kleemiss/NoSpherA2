#pragma once
#include <iosfwd>
#include <string>
#include <vector>

//The literature every method in this program implements, so that a run says which paper its
//numbers come from.  Every entry below was checked against Crossref (api.crossref.org/works/<doi>)
//on 24 Sep 2026: authors, journal, volume, year, first page and DOI all come from the registered
//metadata, not from memory.  The full table, with how each line was verified, is in the vault note
//Software-Notes/NoSpherA2-Codebase/NoSpherA2-Method-Citations-24-Sep-V1.0.
//
//No equation numbers are quoted here on purpose: a wrong one sends a reader to the wrong formula,
//and the papers behind most of these entries are not readable without a subscription.  Where a
//header in this tree does quote one (eli_family.h, Eq. 52/53 of Kohout III), that is the place for it.

namespace citations
{
    //One method per enumerator; a method may carry more than one Reference in the table.
    enum class Method
    {
        NoSpherA2,          //the program itself, and the .tsc format it defined
        HAR,                //Hirshfeld atom refinement / aspherical form factors
        Hirshfeld,          //Hirshfeld (stockholder) partitioning
        IAM,                //Slater/Thakkar independent-atom form factors
        BeckeGrid,          //the multicentre atomic integration grid
        TFVC,               //topological fuzzy Voronoi cell partitioning
        MBIS,               //minimal basis iterative stockholder partitioning
        EMBIS,              //ellipsoidal MBIS - the sigma tensors of make_EMBIS_tensors
        ECP,                //the ECP core correction
        RIFit,              //RI / density-fitted partitioning and form factors
        Embedding,          //electrostatic and multi-layer embedding for HAR
        QTAIM,              //critical points and basin integration
        LIDI,               //localisation and delocalization indices
        ELF,                //electron localization function
        ELID,               //ELI-D
        ELIFamily,          //triplet-coupled ELI-D and ELI-q
        ELIOrbitalFree,     //the PC07 orbital-free ELI estimate
        NCI,                //noncovalent interaction / reduced density gradient plots
        ESP,                //molecular electrostatic potential
        HirshfeldSurface,   //Hirshfeld surfaces, d_i/d_e, fingerprints
        NAONPA,             //natural atomic orbitals and natural population analysis
        NBO,                //the natural bond orbital search
        E2,                 //second-order donor-acceptor perturbation
        NRT,                //natural resonance theory
        NBOProgram,         //the external NBO program and its FILE47 (.47) format
        Molden,             //the molden file format
        SALTED,             //machine-learned electron density
        D4,                 //the D4 dispersion model
        PTB,                //the pTB tight-binding potential
        TRAH,               //trust-region augmented Hessian second-order SCF
        RGBI,               //Roby-Gould bond indices
        Fukui,              //Fukui functions and the dual descriptor
        GordonKim           //Gordon-Kim exchange repulsion
    };

    struct Reference
    {
        Method method;
        const char *tag;   //what the log line is labelled with, e.g. "QTAIM"
        const char *work;  //"Bader, Chem. Rev. 91 (1991) 893"
        const char *doi;   //bare DOI, no https:// prefix
    };

    //Every reference, grouped by method in Method order.
    const std::vector<Reference> &table();

    //"[QTAIM] Bader, Chem. Rev. 91 (1991) 893, DOI 10.1021/cr00005a013"
    std::string format(const Reference &reference);

    //Write one line per reference registered for this method, in the log style of everything else.
    void cite(Method method, std::ostream &os);

    //The same, to std::cout.
    void cite(Method method);
}
