#include "citations.h"

#include <iostream>

namespace citations
{
    const std::vector<Reference> &table()
    {
        //Crossref-confirmed on 24 Sep 2026.  A method with two entries rests on both papers: the
        //method itself and the implementation this code follows.
        static const std::vector<Reference> refs = {
            {Method::NoSpherA2, "NoSpherA2", "Kleemiss et al., Chem. Sci. 12 (2021) 1675", "10.1039/D0SC05526C"},
            {Method::HAR, "HAR", "Jayatilaka & Dittrich, Acta Cryst. A 64 (2008) 383", "10.1107/S0108767308005709"},
            {Method::HAR, "HAR", "Capelli et al., IUCrJ 1 (2014) 361", "10.1107/S2052252514014845"},
            {Method::Hirshfeld, "Hirshfeld", "Hirshfeld, Theor. Chim. Acta 44 (1977) 129", "10.1007/BF00549096"},
            {Method::IAM, "IAM", "Kleemiss, Peyerimhoff & Bodensteiner, J. Appl. Cryst. 57 (2024) 161", "10.1107/S1600576723010981"},
            {Method::IAM, "IAM", "Koga, Tatewaki & Thakkar, Phys. Rev. A 47 (1993) 4510", "10.1103/PhysRevA.47.4510"},
            {Method::BeckeGrid, "Becke grid", "Becke, J. Chem. Phys. 88 (1988) 2547", "10.1063/1.454033"},
            {Method::TFVC, "TFVC", "Salvador & Ramos-Cordoba, J. Chem. Phys. 139 (2013) 071103", "10.1063/1.4818751"},
            {Method::MBIS, "MBIS", "Verstraelen et al., J. Chem. Theory Comput. 12 (2016) 3894", "10.1021/acs.jctc.6b00456"},
            //The E is "ellipsoidal": this is the anisotropic sigma tensor, not a variant name of ours.
            {Method::EMBIS, "EMBIS", "Nielsen & Jensen, J. Chem. Theory Comput. 21 (2025) 8753", "10.1021/acs.jctc.5c00788"},
            {Method::ECP, "ECP", "Kleemiss et al., J. Appl. Cryst. 58 (2025) 374", "10.1107/S1600576725000901"},
            {Method::RIFit, "RI fit", "Seifert et al., Z. Kristallogr. Cryst. Mater. 241 (2026) 283", "10.1515/zkri-2026-0013"},
            {Method::Embedding, "Embedding", "Landeros-Rivera & Kleemiss, J. Appl. Cryst. 59 (2026)", "10.1107/S160057672600717X"},
            {Method::QTAIM, "QTAIM", "Bader, Chem. Rev. 91 (1991) 893", "10.1021/cr00005a013"},
            {Method::LIDI, "LI/DI", "Fradera, Austen & Bader, J. Phys. Chem. A 103 (1999) 304", "10.1021/jp983362q"},
            {Method::ELF, "ELF", "Becke & Edgecombe, J. Chem. Phys. 92 (1990) 5397", "10.1063/1.458517"},
            {Method::ELID, "ELI-D", "Kohout, Pernal, Wagner & Grin, Theor. Chem. Acc. 112 (2004) 453", "10.1007/s00214-004-0615-y"},
            {Method::ELIFamily, "ELI family", "Kohout, Wagner & Grin, Theor. Chem. Acc. 119 (2008) 413", "10.1007/s00214-007-0396-1"},
            {Method::ELIFamily, "ELI family", "Kohout, Faraday Discuss. 135 (2007) 43", "10.1039/b605951c"},
            {Method::ELIOrbitalFree, "ELI (PC07)", "Perdew & Constantin, Phys. Rev. B 75 (2007) 155109", "10.1103/PhysRevB.75.155109"},
            {Method::NCI, "NCI", "Johnson et al., J. Am. Chem. Soc. 132 (2010) 6498", "10.1021/ja100936w"},
            {Method::ESP, "ESP", "Bonaccorsi, Scrocco & Tomasi, J. Chem. Phys. 52 (1970) 5270", "10.1063/1.1672775"},
            {Method::HirshfeldSurface, "Hirshfeld surface", "Spackman & Jayatilaka, CrystEngComm 11 (2009) 19", "10.1039/B818330A"},
            {Method::HirshfeldSurface, "Hirshfeld surface", "McKinnon, Spackman & Mitchell, Acta Cryst. B 60 (2004) 627", "10.1107/S0108768104020300"},
            {Method::NAONPA, "NAO/NPA", "Reed, Weinstock & Weinhold, J. Chem. Phys. 83 (1985) 735", "10.1063/1.449486"},
            {Method::NBO, "NBO", "Foster & Weinhold, J. Am. Chem. Soc. 102 (1980) 7211", "10.1021/ja00544a007"},
            {Method::E2, "E2", "Reed, Curtiss & Weinhold, Chem. Rev. 88 (1988) 899", "10.1021/cr00088a005"},
            {Method::NRT, "NRT", "Glendening & Weinhold, J. Comput. Chem. 19 (1998) 593", "10.1002/(SICI)1096-987X(19980430)19:6<593::AID-JCC3>3.0.CO;2-M"},
            {Method::NRT, "NRT", "Glendening & Weinhold, J. Comput. Chem. 19 (1998) 610", "10.1002/(SICI)1096-987X(19980430)19:6<610::AID-JCC4>3.0.CO;2-U"},
            {Method::NBOProgram, "NBO 7", "Glendening, Landis & Weinhold, J. Comput. Chem. 40 (2019) 2234", "10.1002/jcc.25873"},
            {Method::Molden, "molden", "Schaftenaar & Noordik, J. Comput.-Aided Mol. Des. 14 (2000) 123", "10.1023/A:1008193805436"},
            {Method::SALTED, "SALTED", "Grisafi, Lewis, Rossi & Ceriotti, J. Chem. Theory Comput. 19 (2023) 4451", "10.1021/acs.jctc.2c00850"},
            {Method::D4, "D4", "Caldeweyher et al., J. Chem. Phys. 150 (2019) 154122", "10.1063/1.5090222"},
            {Method::PTB, "pTB", "Grimme, Muller & Hansen, J. Chem. Phys. 158 (2023) 124111", "10.1063/5.0137838"},
            {Method::TRAH, "TRAH", "Helmich-Paris, J. Chem. Phys. 154 (2021) 164104", "10.1063/5.0040798"},
            {Method::RGBI, "RGBI", "Roby, Mol. Phys. 27 (1974) 81", "10.1080/00268977400100071"},
            {Method::RGBI, "RGBI", "Gould et al., Theor. Chem. Acc. 119 (2008) 275", "10.1007/s00214-007-0282-x"},
            {Method::Fukui, "Fukui", "Parr & Yang, J. Am. Chem. Soc. 106 (1984) 4049", "10.1021/ja00326a036"},
            {Method::Fukui, "Fukui", "Morell, Grand & Toro-Labbe, J. Phys. Chem. A 109 (2005) 205", "10.1021/jp046577a"},
            {Method::GordonKim, "Gordon-Kim", "Gordon & Kim, J. Chem. Phys. 56 (1972) 3122", "10.1063/1.1677649"},
        };
        return refs;
    }

    std::string format(const Reference &reference)
    {
        return "[" + std::string(reference.tag) + "] " + reference.work + ", DOI " + reference.doi;
    }

    void cite(Method method, std::ostream &os)
    {
        for (const Reference &reference : table())
            if (reference.method == method)
                os << format(reference) << std::endl;
    }

    void cite(Method method)
    {
        cite(method, std::cout);
    }

    namespace
    {
        //Readers run one after another on one thread, so a plain vector is enough.
        std::vector<Method> &queued()
        {
            static std::vector<Method> pending;
            return pending;
        }
    }

    void queue(Method method)
    {
        queued().push_back(method);
    }

    void flush(std::ostream &os)
    {
        for (const Method method : queued())
            cite(method, os);
        queued().clear();
    }
}
