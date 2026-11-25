#include "pch.h"
#include "convenience.h"
#include "git_commit.h"

std::string NoSpherA2_message(bool no_date)
{
	std::string t = "    _   __     _____       __              ___   ___\n";
	t.append("   / | / /___ / ___/____  / /_  ___  _____/   | |__ \\\n");
	t.append("  /  |/ / __ \\\\__ \\/ __ \\/ __ \\/ _ \\/ ___/ /| | __/ /\n");
	t.append(" / /|  / /_/ /__/ / /_/ / / / /  __/ /  / ___ |/ __/\n");
	t.append("/_/ |_/\\____/____/ .___/_/ /_/\\___/_/  /_/  |_/____/\n");
	t.append("                /_/\n");
	t.append("This software is part of the cuQCT software suite developed by Florian Kleemiss.\n");
	t.append("Please give credit and cite corresponding pieces!\nThis Software is published with BSD-2 clause license.\n");
	if (!no_date)
	{
		t.append("List of contributors of pieces of code or functionality:\n");
		t.append("      Florian Kleemiss,\n");
		t.append("      Emmanuel Hupf,\n");
		t.append("      Alessandro Genoni,\n");
		t.append("      Lukas M. Seifert,\n");
		t.append("      Daniel Bruex,\n");
		t.append("      Marti Gimferrer,\n");
		t.append("      Anker Nielsen,\n");
		t.append("      Lucas Militao,\n");
		t.append("      Johannes Bartusel,\n");
		t.append("      and many more in communications or by feedback!\n");
		t.append("NoSpherA2 uses featomic, Metatensor, and the mdspan library, as well as OCC for the calculation of wavefunctions, when required.\n");
		t.append("The used packages are published under BSD-3 clause License or explicit consent for the use in this project was given.\n");
		t.append("Please see, respectively:\n");
		t.append("   https://github.com/Luthaf/featomic\n");
		t.append("   https://github.com/lab-cosmo/metatensor\n");
#if defined(NSA2_ARMPL)
		t.append("   This software utilizes Arm Performance Libraries (https://developer.arm.com/Tools%20and%20Software/Arm%20Performance%20Libraries) for optimized mathematical computations\n");
#elif defined(NSA2_OPENBLAS)
		t.append("   This software utilizes OpenBLAS (https://github.com/OpenMathLib/OpenBLAS) for optimized mathematical computations\n");
#else
		t.append("   This software utilizes Intel(c) Math Kernel Library (oneMKL), version 2025.2.0.629, for optimized mathematical computations\n");
#endif
		t.append("OCC can be found at: https://github.com/peterspackman/occ\n");
		t.append("OCC was published at        : Spackman, J. Open Source Softw. 2026, 11, 9609, doi:10.21105/joss.09609.\n");
		t.append("NoSpherA2 was published at  : Kleemiss et al. Chem. Sci., 2021, 12, 1675 - 1692.\n");
		t.append("Slater IAM was published at : Kleemiss et al. J. Appl. Cryst. 2024, 57, 161 - 174.\n");
		t.append("ECP correction functions at : Kleemiss et al. J. Appl. Cryst. 2025, 58, 374 - 382.\n");
		t.append("Aux basis /RI partitioning  : Seifert et al. Z. Krist. - Cryst. Mat. 2026, 241, 283 - 295.\n");
		t.append("Embedding for HAR at        : Landeros-Rivera & Kleemiss, J. Appl. Cryst. 2026, 59, 10.1107/S160057672600717X.\n");
		t.append("TFVC partitioning at        : Gimferrer et al. TBA.\n");
		t.append("MBIS partitioning at        : Verstraelen et al. J. Chem. Theory Comput. 2016, 12, 3894 - 3912.\n");
		t.append("EMBIS partitioning at       : Nielsen & Jensen, J. Chem. Theory Comput. 2025, 21, 8753 - 8761.\n");
	}
	return t;
}

std::string WAEmpiRe_message()
{
	std::string t = " _   _   __ __    _____               __ ___  \n";
	t.append("| | / | / /   |  / ___/_   __  ____  /_/ __ \\___ \n");
	t.append("| |/  |/ / /| | / /__/  | /  |/ __ \\__/ /_/ / _ \\ \n");
	t.append("|  / |  / __  |/ /__/ /||//| / /_/ / / __  /  __/ \n");
	t.append("|_/  |_/_/  |_/____/_/ |_/ |/ .___/_/_/  |_\\___/ \n");
	t.append("                           /_/				\n");
	t.append("Wavefunction Augmentation through Emprirical Restraints \n");
	t.append("is part of the NoSpherA2 project developed by Florian Kleemiss. \n");
	return t;
}



std::string build_date = ("This Executable was built on: " + std::string(__DATE__) + " " + std::string(__TIME__) + " from git commit " NOSPHERA2_GIT_COMMIT "\n");
