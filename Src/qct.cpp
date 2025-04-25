#include "pch.h"
#include "convenience.h"
#include "fchk.h"
#include "basis_set.h"
#include "wfn_class.h"
#include "b2c.h"
#include "bondwise_analysis.h"

int QCT(options& opt)
{
	using namespace std;
    bool expert = false;
	bool end = false;
	vector<WFN> wavy;
	int activewave=0;
	int cpus = -1;
	string wavename;
	int ncpus=0;
	float mem=0.0;
	unsigned int counter = 0;
	cls();
	while(!end){
		char sel;
		std::cout << "This Executable was built on: " + string(__DATE__) + " " + string(__TIME__) + "\n";
		std::cout << "What do you want to do?\n";
		if(wavy.size() > 0) {
			std::cout << "Active Wavefunction: " << wavy[activewave].get_path();
			string temp = wavy[activewave].get_basis_set_name();
			if(temp.length()>2)std::cout << " (" << temp << ") ";
			else std::cout << " (no basis set loaded)";
			if(opt.debug)std::cout << temp;
			if(wavy[activewave].get_modified())std::cout << "*";
			if(expert)std::cout << " EXPERT MODE!";
			std::cout << endl;
		}
		std::cout << "R) Read in a new file" << endl;
		if(wavy.size() > 0){
			std::cout << "M) Modify an already loaded wavefunction"<<endl 
				<< "S) Save the active wavefunction";
			if (wavy[activewave].cub.size() > 0)
				std::cout << "or cube(s)" << endl;
			else std::cout << endl;
			std::cout << "O) Sort the exponents in the wavefunction" << endl
				<< "N) Start NCI/Cube calculation plugin" << endl
				<< "B) Read a basis set" << endl
				<< "U) Check the unit of the atom positions" << endl;
			if(wavy.size() > 1)std::cout << "A) Activate another wavefunction" << endl;
			else std::cout << "-) Activate another wavefunction" << endl;
			if(wavy[activewave].cub.size()>0)std::cout << "C) Work with cube files loaded" << endl;
			else std::cout << "-) Work with cube files loaded" << endl;
		}
		else{
			std::cout << "-) Modify an already loaded wavefunction"<<endl
				<< "-) Save the active wavefunction" << endl
				<< "-) Sort the exponents in the wavefunction" << endl
				<< "-) Start NCI/Cube calculation plugin" << endl
				<< "-) Read a basis set" << endl
				<< "-) Check the unit of the atom positions" << endl;
			if(wavy.size() > 1)std::cout << "A) Activate another wavefunction" << endl;
			else std::cout << "-) Activate another wavefunction" << endl;
			if(wavy.size() > 0 && wavy[activewave].cub.size()>0)std::cout << "C) Work with cube files loaded" << endl;
			else std::cout << "-) Work with cube files loaded" << endl;
		}
		std::cout << "E) Toggle Expert mode (Disable assumptions)" << endl
			 << "L) Limit the number of CPUs being used. Current value (-1 corresponds to all): " << cpus << endl
			 << "Q) Quit the program" << endl;
		std::cin >> sel;
		cls();
		vector < vector <unsigned int> > selection;
		switch(sel){
			case 'c':
			case 'C':{
				unsigned int nr_cubes=0;
				selection.resize(2);
				for(int w=0; w<wavy.size(); w++) for(int c=0; c<wavy[w].cub.size(); c++) nr_cubes++;
				std::cout << "What do you want to do?" << endl
					<< "1) Perform mathematic operation on one cube" << endl;
				if(nr_cubes>1)std::cout << "2) Perform mathematic operation on two cubes" << endl;
				std::cout << "3) Analyze a cube" << endl
					<< "0) Go back to main menu" << endl;
				int temp;
				std::cin >> temp;
				cls();
				switch(temp){
					case 1:{
						if(nr_cubes<=0){
							std::cout << "No cubes loaded!" << endl;
							break;
						}
						std::cout << "Which operation do you want to perform?" << endl
							 << "1) Integrate values inside a cube" << endl
							 << "2) Periodically replicate a cube" << endl
							 << "3) Apply a threshhold to a cube" << endl
							 << "0) Get back to previous menu" << endl;
						int sel=0;
						std::cin >> sel;
						selection[0].resize(1);
						selection[1].resize(1);
						if(nr_cubes>=1&&sel!=0)	select_cubes(selection, wavy, 1);
						if(wavy[selection[0][0]].cub[selection[1][0]].get_loaded()==false){
							if (opt.debug)std::cout << "Loading full file now!" << endl;
							wavy[selection[0][0]].cub[selection[1][0]].read_file(true,false,wavy[activewave]);
						}
						switch(sel){
						case 1:{
							cls();
							std::cout << "Integrated value of this cube is: ";
							std::cout << scientific << setprecision(8) << setw(16) << wavy[selection[0][0]].cub[selection[1][0]].sum() << endl;
							break;
						}
						case 2:{
							string result=wavy[selection[0][0]].cub[selection[1][0]].super_cube(wavy[selection[0][0]]);
							cls();
							std::cout << "New cube saved as: " << result << endl;
							break;
						}
						case 3:{
							double thresh=0.0;
							std::cout << "Please enter threshhold: ";
							std::cin >> thresh;
							wavy[selection[0][0]].cub[selection[1][0]].thresh(thresh);
							break;
						}
						case 4:{
								cls();
								std::cout << "Integrated absolute_differnce value of this cube is: ";
								std::cout << scientific << setprecision(8) << setw(16) << wavy[selection[0][0]].cub[selection[1][0]].diff_sum() << endl;
								break;
						}
                        case 0:{
                            end=true;
                            break;
                        }
						default:
                           std::cout << "Sorry, i did not get that ..." << endl;
							break;
						}
						break;
					}
					case 2:{
						if(nr_cubes>1){
							while(!end){
								std::cout << "Which operation do you want to perform? ";
								if(expert)std::cout << "(The 1. cube selected will be altered by the 2. one in cases of x=)";
								std::cout << endl
									<< "1) +" << endl
									<< "2) -" << endl
									<< "3) *" << endl
									<< "4) /" << endl
									<< "5) RSR (Real space R-Value)" << endl
									<< "6) Mask (if(value!=0) value, else 0)" << endl
									<< "7) Negative mask (if(value==0) value, else 0" << endl
									<< "8) Mask with threshhold (if value(2) < thresh = 0, else = value(1))" << endl
									<< "9) Weighted Jaccord similarity and distance" << endl;
								if(expert)std::cout << "11) +=" << endl
										<< "12) -=" << endl
										<< "13) *=" << endl
										<< "14) /=" << endl;
								std::cout << "0) Go back to previous menu" << endl;
								std::cin >> temp;
								selection[0].resize(2);
								selection[1].resize(2);
								if(nr_cubes>=2&&temp!=0){
									select_cubes(selection, wavy, 2);
									if(opt.debug){
										std::cout << "Selection:" << endl;
										for(int i=0; i<2; i++)std::cout << selection[0][i] << "." << selection[1][i] << endl;
									}
								}
								cls();
								switch(temp){
									case 1:{
										wavy[selection[0][0]].push_back_cube(wavy[selection[0][0]].cub[selection[1][0]]+wavy[selection[0][1]].cub[selection[1][1]]);
										if(wavy[selection[0][0]].cub[wavy[selection[0][0]].cub.size()-1].get_size(0)!=0)std::cout << "Operation succesfull!" << endl;
										else std::cout << "Sorry, something went wrong!" << endl;
									}
										break;
									case 2:{
										wavy[selection[0][0]].push_back_cube(wavy[selection[0][0]].cub[selection[1][0]]-wavy[selection[0][1]].cub[selection[1][1]]);
										if(wavy[selection[0][0]].cub[wavy[selection[0][0]].cub.size()-1].get_size(0)!=0)std::cout << "Operation succesfull!" << endl;
										else std::cout << "Sorry, something went wrong!" << endl;
									}
										break;
									case 3:{
										wavy[selection[0][0]].push_back_cube(wavy[selection[0][0]].cub[selection[1][0]]*wavy[selection[0][1]].cub[selection[1][1]]);
										if(wavy[selection[0][0]].cub[wavy[selection[0][0]].cub.size()-1].get_size(0)!=0)std::cout << "Operation succesfull!" << endl;
										else std::cout << "Sorry, something went wrong!" << endl;
									}
										break;
									case 4:{
										wavy[selection[0][0]].push_back_cube(wavy[selection[0][0]].cub[selection[1][0]]/wavy[selection[0][1]].cub[selection[1][1]]);
										if(wavy[selection[0][0]].cub[wavy[selection[0][0]].cub.size()-1].get_size(0)!=0)std::cout << "Operation succesfull!" << endl;
										else std::cout << "Sorry, something went wrong!" << endl;
									}
										break;
									case 5:{
										for(int i=0; i<2; i++)
											if(wavy[selection[0][i]].cub[selection[1][i]].get_loaded()==false){
												if (opt.debug)std::cout << "Loading full file now!" << endl;
												wavy[selection[0][i]].cub[selection[1][i]].read_file(true,false,wavy[activewave]);
											}
										double result=wavy[selection[0][0]].cub[selection[1][0]].rrs(wavy[selection[0][1]].cub[selection[1][1]]);
										if(result!=-1){
											std::cout << "Operation succesfull!" << endl;
											std::cout << "RSR: " << result << endl;
										}
										else std::cout << "Sorry, something went wrong!" << endl;
										}
										break;
									case 6: {
										for (int i = 0; i < 2; i++)
											if (wavy[selection[0][i]].cub[selection[1][i]].get_loaded() == false) {
												if (opt.debug)std::cout << "Loading full file now!" << endl;
												wavy[selection[0][i]].cub[selection[1][i]].read_file(true, false, wavy[activewave]);
											}
										if (wavy[selection[0][0]].cub[selection[1][0]].mask(wavy[selection[0][1]].cub[selection[1][1]]))std::cout << "Operation succesfull!" << endl;
										else std::cout << "Sorry, something went wrong!" << endl;
										}
										break;
									case 7:
										for(int i=0; i<2; i++)
											if(wavy[selection[0][i]].cub[selection[1][i]].get_loaded()==false){
												if (opt.debug) 
													std::cout << "Loading full file now!" << endl;
												wavy[selection[0][i]].cub[selection[1][i]].read_file(true,false,wavy[activewave]);
											}
										if(wavy[selection[0][0]].cub[selection[1][0]].negative_mask(wavy[selection[0][1]].cub[selection[1][1]])) 
											std::cout << "Operation succesfull!" << endl;
										else 
											std::cout << "Sorry, something went wrong!" << endl;
										break;
									case 8:{
										for(int i=0; i<2; i++)
											if(wavy[selection[0][i]].cub[selection[1][i]].get_loaded()==false){
												if (opt.debug) 
													std::cout << "Loading full file now!" << endl;
												wavy[selection[0][i]].cub[selection[1][i]].read_file(true,false,wavy[activewave]);
											}
										std::cout << "Please give threshhold to use: ";
										double thresh=0.0;
										std::cin >> thresh;
										if(wavy[selection[0][0]].cub[selection[1][0]].mask(wavy[selection[0][1]].cub[selection[1][1]],thresh)) 
											std::cout << "Operation succesfull!" << endl;
										else 
											std::cout << "Sorry, something went wrong!" << endl;
										break;
									}
									case 9:{
										for(int i=0; i<2; i++)
											if(wavy[selection[0][i]].cub[selection[1][i]].get_loaded()==false){
												if (opt.debug)std::cout << "Loading full file now!" << endl;
												wavy[selection[0][i]].cub[selection[1][i]].read_file(true,false,wavy[activewave]);
											}
										double result=wavy[selection[0][0]].cub[selection[1][0]].jaccard(wavy[selection[0][1]].cub[selection[1][1]]);
										if(result!=-1){
											std::cout << "Operation succesfull!" << endl;
											std::cout << "Jaccord similarity: " << result << " Jaccord distance: " << 1-result << endl;
										}
										else 
											std::cout << "Sorry, something went wrong!" << endl;
									}
									break;
									case 11:
										for(int i=0; i<2; i++)
											if(wavy[selection[0][i]].cub[selection[1][i]].get_loaded()==false){
												if (opt.debug)std::cout << "Loading full file now!" << endl;
												wavy[selection[0][i]].cub[selection[1][i]].read_file(true,false,wavy[activewave]);
											}
										if(wavy[selection[0][0]].cub[selection[1][0]]+=wavy[selection[0][1]].cub[selection[1][1]]) 
											std::cout << "Operation succesfull!" << endl;
										else 
											std::cout << "Sorry, something went wrong!" << endl;
										break;
									case 12:
										for(int i=0; i<2; i++)
											if(wavy[selection[0][i]].cub[selection[1][i]].get_loaded()==false){
												if (opt.debug)std::cout << "Loading full file now!" << endl;
												wavy[selection[0][i]].cub[selection[1][i]].read_file(true,false,wavy[activewave]);
											}
										if(wavy[selection[0][0]].cub[selection[1][0]]-=wavy[selection[0][1]].cub[selection[1][1]])std::cout << "Operation succesfull!" << endl;
										else std::cout << "Sorry, something went wrong!" << endl;
										break;
									case 13:
										for(int i=0; i<2; i++)
											if(wavy[selection[0][i]].cub[selection[1][i]].get_loaded()==false){
												if (opt.debug)std::cout << "Loading full file now!" << endl;
												wavy[selection[0][i]].cub[selection[1][i]].read_file(true,false,wavy[activewave]);
											}
										if(wavy[selection[0][0]].cub[selection[1][0]]*=wavy[selection[0][1]].cub[selection[1][1]])std::cout << "Operation succesfull!" << endl;
										else std::cout << "Sorry, something went wrong!" << endl;
										break;
									case 14:
										for(int i=0; i<2; i++)
											if(wavy[selection[0][i]].cub[selection[1][i]].get_loaded()==false){
												if (opt.debug)std::cout << "Loading full file now!" << endl;
												wavy[selection[0][i]].cub[selection[1][i]].read_file(true,false,wavy[activewave]);
											}
										if(wavy[selection[0][0]].cub[selection[1][0]]/=wavy[selection[0][1]].cub[selection[1][1]])std::cout << "Operation succesfull!" << endl;
										else std::cout << "Sorry, something went wrong!" << endl;
										break;
									case 0:
										end=true;
										break;
									default:
										std::cout << "Sorry, i did not get that..." << endl;
										break;
								}
							}
							end=false;
						}
						break;
					}
					case 3:
						std::cout << "Which analysis do you want to perform?" << endl
							<< "1) Separate into Basins according to critical point search" << endl;
							//and try new BCP implementation"
						if(expert)std::cout << "2) Separate into Basins according to critical point search and try new BCP implementation" << endl;
						std::cout << "0) Go back to selection" << endl;
						std::cin >> temp;
						selection[0].resize(1);
						selection[1].resize(1);
						if(nr_cubes>=1&&temp!=0) select_cubes(selection, wavy, 1, false, opt.debug);
						if(opt.debug) std::cout << "selection: " << selection[0][0] << " " << selection[1][0] << endl;
						switch(temp){
							case 1:
								if(wavy[selection[0][0]].cub[selection[1][0]].get_loaded()==false){
									std::cout << "Loading full file now!" << endl;
									if(wavy[selection[0][0]].cub[selection[1][0]].read_file(true,false,wavy[activewave])==false){
										std::cout << "ERROR reading full file! Aborting" << endl;
										break;
									}
								}
								if(b2c(wavy[selection[0][0]].cub[selection[1][0]], wavy[selection[0][0]].atoms, opt.debug,false)==false)std::cout << "something went wrong!" << endl;
								break;
							case 2:
								if(expert){
									if(wavy[selection[0][0]].cub[selection[1][0]].get_loaded()==false){
										std::cout << "loading full file now!" << endl;
										if(wavy[selection[0][0]].cub[selection[1][0]].read_file(true,false,wavy[activewave])==false){
											std::cout << "ERROR reading full file! Aborting!" << endl;
											break;
										}
									}
									if(b2c(wavy[selection[0][0]].cub[selection[1][0]], wavy[selection[0][0]].atoms, opt.debug,true)==false)std::cout << "something went wrong" << endl;
									break;
								}
								break;
							case 0:
								break;
							default:
								std::cout << "sorry, i did not understand that, please try again!" << endl;
								break;
						}
					}
					break;
				}		
			case 'r':
			case 'R':{
				if(expert){
					std::cout << "What kind of file do you want to load?" << endl
							<< "W) Wavefunction in known format (gbw/wfn/wfx/fchk/molden)"<< endl
							<< "G) Grid or Cube file" << endl
							<< "E) Exit to main menu" << endl;
					std::cin >> sel;
					cls();
					switch(sel){
						case 'G':
						case 'g':{
							std::filesystem::path path;
							while (!exists(path)) {
								std::cout << "Please give Path to the cube you want to read: ";
								std::cin >> path;
								if (exists(path)) break;
								else std::cout << "Sorry, couldn't find the file! Try again!" << endl;
							}
							if(wavy.size() == 0){
								activewave=0;
								wavy.push_back(WFN(path));
							}
							wavy[activewave].push_back_cube(cube(path, false, wavy[activewave], std::cout, expert);
							cls();
							break;
						}
						case 'W':
						case 'w':{
							while(!end){
								string path;
								bool new_wfn=false;
								std::cout << "Path to the wavefunction file: ";
								std::cin >> path;
								int tries=0;
								if(path.find(".wfn")==-1&&path.find(".ffn")==-1){
									do {
										std::cout << "This doesn't look like a .wfn or .ffn file, try again: " << endl;
										tries++;
										std::cin >> path;
									} while (path.find(".wfn")==-1&&path.find(".ffn")==-1&&tries<3);
									if(tries == 3) {
										std::cout << "Sorry, make sure you know the filename and try again!" << endl;
										break;
									}
								}
								if(activewave<wavy.size()&&wavy.size()>0){
									std::cout << "This will delete the previously loaded wavefunction "<< wavy[activewave].get_path() <<"! Are you sure?";
									if(!yesno()){
										std::cout << "Do you want to load it as a new wavefunction? " << endl;
										if(yesno()) new_wfn=true;
										end=true;
										break;
									}
									wavy.erase(wavy.begin()+activewave);
									if(path.find(".wfn")!=-1) wavy.insert(wavy.begin()+activewave-1,WFN(2));
									else wavy.insert(wavy.begin()+activewave-1,WFN(4));
								}
								else if(wavy.size()==0||activewave==wavy.size()||new_wfn) {
									if(opt.debug) std::cout << "Making a new wavefunction!" << endl;
									if(path.find(".wfn")!=-1) wavy.push_back(WFN(2));
									else wavy.push_back(WFN(4));
									if(wavy.size()>1) activewave++;
									if(opt.debug) std::cout << "Size: " << wavy.size() << " active: " << activewave << endl;
								}
								wavy[activewave].read_known_wavefunction_format(path, std::cout, opt.debug);
								cls();
								if(opt.debug) std::cout << wavy[activewave].get_ncen () << endl;
								break;
							}
							cls();
							end=false;
							break;
						}
						case 'e':
						case 'E':
							break;
						default:
							std::cout << "Sorry, did not understand that, going back to main menu!" << endl;
							break;
					}
				}
				else{
					string filename;
					vector <string> temp;
					temp.resize(4);
					temp[0]="Wavefunction files (wfn,ffn,fchk,wfx) | *.wfn *.ffn *.Fchk *.fchk *.FChk *.wfx";
					temp[1]="Cube files (cub, cube, grd) | *.cube *.cub *.grd";
					temp[2]="Crystal Output (out) | *.out";
					temp[3]="All filetypes | *";
					if(!open_file_dialog(filename,opt.debug,temp)){
						std::cout << "Error encountered!" << endl;
						break;
					}
					else{
                        wavy[activewave].read_known_wavefunction_format(filename, std::cout, opt.debug);
					}
				}
				break;
			}
			case 'M':
			case 'm':{
				if(wavy.size()<1){
					std::cout << "First you need to read a wavefunction!" << endl;
					break;
				}
				int msel=0;
				string label="?";
				float x,y,z=0.0;
				int charge=0;
				std::cout << "What do you want to do?"<<endl
					<< "D) Delete a center or set of values for a certain function"<<endl
					<< "A) Add a center" << endl
					<< "C) Change a value" << endl
					<< "E) Exit to Main Menu" << endl;
				std::cin >> sel;
				cls();
				switch(sel){
					case 'D':
					case 'd':
						std::cout << "A) Atom"<< endl
							<< "E) Exponent" << endl
							<< "B) Back to main menu"<<endl;
						std::cin >> sel;
						end=false;
						switch(sel){
							case 'A':
							case 'a':
								std::cout << "The list of centers:\n";
								wavy[activewave].list_centers();
								std::cout << "Which center do you want to delete?\n";
								std::cin >> msel;
								if (msel>=0){
									if(wavy[activewave].remove_center(msel)){
										std::cout << "Deleted center nr " << msel << " succesfully! Going back to Main Menu.\n";
										wavy[activewave].set_modified ();
									}
									else std::cout << "Something went wrong... Sorry, let's start again...\n";
								}
								else std::cout << "Wrong selection, start again!\n";
								break;
							case 'E':
							case 'e':
								while(!end){
									wavy[activewave].list_primitives ();
									std::cout << "Which exponent out of "<< wavy[activewave].get_nex() << " do you want to delete?\n";
									std::cin >> msel;
									if(msel<0||msel>wavy[activewave].get_nex()){
										std::cout << "Sorry, wrong input";
										continue;
									}
									std::cout << "This is the set of information you requested:\n";
									wavy[activewave].print_primitive(msel);
									end=true;
								}
								break;
							case 'B':
							case 'b':
								cls();
								break;
						}
						break;
					case 'A':
					case 'a':
						std::cout << "A) Atom\nE) Exponent\nB) Back to main menu" << endl;
						std::cin >> sel;
						end=false;
						switch(sel){
							case 'E':
							case 'e':
								while(!end){
									std::cout << "Please remember that this exponent will be apended to the data structure, "
										<< "it will not be sorted in any way!\nCentre Assignement: ";
									int temp_cen=0;
									std::cin >> temp_cen;
									if (temp_cen>wavy[activewave].get_ncen ()||temp_cen<0){
										std::cout << "Wrong input, start over!";
										continue;
									}
									std::cout << "Type: ";
									int temp_type=0;
									std::cin >> temp_type;
									std::cout << "Exponent: ";
									double temp_exp=0.0;
									std::cin >> temp_exp;
									vector<double> temp_val;
									temp_val.resize(wavy[activewave].get_nmo());
									for (int i=0; i<wavy[activewave].get_nmo();i++){
										std::cout << "Enter coefficient for MO " << i << ":";
										std::cin >> temp_val[i];
										if(temp_val[i]<-1000||temp_val[i]>1000){
											std::cout << "Wrong input, please try again...\n";
											i--;
											continue;
										}
									}   
									std::cout << "Let me recapitulate: Center " << temp_cen << " type: " << temp_type << " exp: " << temp_exp
										<< " and the MO coefficients:\n";
									for(int i=0; i<wavy[activewave].get_nmo(); i++){
										std::cout << temp_val[i] <<"   ";
										if( i%5==0)std::cout << endl;
									}
									std::cout << "is this correct?";
									if(yesno()) end=true;
								}
								break;
							case 'A':
							case 'a':
								while(!end){
									std::cout << "The list of centers:\n";
									wavy[activewave].list_centers();
									std::cout << "Which center do you want to add?\nlabel:";
									std::cin >> label;
									std::cout << "x: ";
									std::cin >> x;
									if(x < -99.999 || x > 99.999){
										std::cout << "Sorry, number too large\n";
										continue;
									}   
									std::cout << "y: ";
									std::cin >> y;
									if(y < -99.999 || y > 99.999){
										std::cout << "Sorry, number too large\n";
										continue;
									}
									std::cout << "z: ";
									std::cin >> z;
									if(z < -99.999 || z > 99.999){
										std::cout << "Sorry, number too large\n";
										continue;
									}
									std::cout << "charge: ";
									std::cin >> charge;
									if(charge <= 0 || charge > 118){
										std::cout << "Sorry, that atom is not yet disovered\n";
										continue;
									}
									wavy[activewave].push_back_atom(label, x, y, z, charge);
								}
								break;
							case 'B':
							case 'b':
								cls();
								break;
						}
						break;
					case 'C':
					case 'c':
						if(wavy.size() <1) continue;
						std::cout << "The center/type/exponent status until now is:\n";
						wavy[activewave].list_primitives ();
						std::cout << "What do you want to change?\nC) Center Assignement\nT) Type assignement\nE) Exponent\nM) MO coefficient\nQ) Quit\n";
						std::cin >> sel;
						switch (sel){
							case 'Q':
							case 'q':
								cls();
								break;
							case 'C':
							case 'c':
								std::cout << "Which one do you want to change? (0=return to menu)";
								std::cin >> msel;
								if(msel>wavy[activewave].get_nex()||msel<0){ 
									std::cout << "Wrong input, start again!\n";
									break;
								}
								else if(msel==0) break;
								wavy[activewave].change_center (msel);
								break;
							case 'E':
							case 'e':
								std::cout << "What is the nr. of the exponent you want to change? (0=return to menu)" << endl;
								std::cin >> msel;
								if(msel>wavy[activewave].get_nex()||msel<0){ 
									std::cout << "Wrong input, start again!\n";
									break;
								}
								else if(msel==0) break;
								wavy[activewave].change_exponent (msel);
								break;
							case 'T':
							case 't':
								std::cout << "What is the nr. of the type you want to change? (0=return to menu)" << endl;
								std::cin >> msel;
								if(msel>wavy[activewave].get_nex()||msel<0){ 
									std::cout << "Wrong input, start again!\n";
									break;
								}
								else if(msel==0) break;
								wavy[activewave].change_type (msel);
								break;
							case 'M':
							case 'm':
								bool end=false;
								while(!end){
									std::cout << "Which MO out of " << wavy[activewave].get_nmo() << " MOs? (0=return to menu)\n";
									int MOsel;
									std::cin >> MOsel;
									if(MOsel>wavy[activewave].get_nmo()||MOsel<0){
										std::cout << "This is not a valid choice...\n";
										continue;
									}
									else if(MOsel==0) break;
									else{
										wavy[activewave].change_MO_coef(MOsel);
										Enter();
										end=true;
									}
								}
							break;
						}
					case 'E':
					case 'e':
						break;
				}
				break;
			}
			case 'n':
			case 'N':{
				std::cout << "Bondwise analysis (B) or NCIplot features (N)? ";
				char seln;
				std::cin >> seln;
				switch(seln){
				case 'B':
				case 'b':{
					string inputfile;
					std::cout << "Which file to use for definition of bonds?" << endl;
					std::cin >> inputfile;
					if(autobonds(opt.debug, wavy[activewave],inputfile,false)!=1) std::cout << "Sorry, looks like something went wrong..." << endl;
					break;
				}
				case 'N':
				case 'n':{
					std::cout << "Starting acuNCI..." << endl;
					acu_nci(wavy,opt.debug,expert);
					break;
				}
				default:
					std::cout << "sorry, didn't understand that" << endl;
					break;
				}
				//cls();
				break;
			}
			case 'S':
			case 's': {
				if (wavy.size() < 1) {
					cls();
					continue;
				}
				vector <string> endings;
				bool w_or_c = true;
				bool convert = false;
				if (wavy[activewave].cub.size() > 0 && wavy[activewave].get_origin() != 3) {
					while (true) {
						std::cout << "Do you want to save the wavefunction (W) or associated cubes (C) or convert cubes into non-cube format (N)? ";
						char input;
						std::cin >> input;
						switch (input) {
						case 'C':
						case 'c': {
							w_or_c = false;
							break;
						}
						case 'W':
						case 'w': {
							w_or_c = true;
							break;
						}
						case 'n':
						case 'N': {
							convert = true;
							break;
						}
						default:
							std::cout << "Sorry, i did not get that! Try again!" << endl;
							continue;
						}
						break;
					}
				}
				else if (wavy[activewave].get_origin() == 3) {
					while (true) {
						std::cout << "Do you want to save the cube in .cube (C) format or convert cubes into non-cube format (N)? ";
						char input;
						std::cin >> input;
						switch (input) {
						case 'C':
						case 'c': {
							convert = false;
							break;
						}
						case 'n':
						case 'N': {
							convert = true;
							break;
						}
						default:
							std::cout << "Sorry, i did not get that! Try again!" << endl;
							continue;
						}
						break;
					}
				}

				if (convert) {
					while (true) {
						std::cout << "Which format do you want to convert to?\nD) DGrid\nX) XD-Graph\n";
						char input;
						std::cin >> input;
						switch (input) {
						case 'D':
						case 'd': {
							if (wavy[activewave].cub.size() >= 1) {
								std::cout << "Which cube do you want to save?" << endl;
								selection.resize(2);
								selection[0].resize(1);
								selection[1].resize(1);
								select_cubes(selection, wavy, 1);
							}
							string path = wavy[selection[0][0]].cub[selection[1][0]].get_path();
							path = path.substr(0, path.find(".cub"));
							if (!wavy[selection[0][0]].cub[selection[1][0]].get_loaded())
								wavy[selection[0][0]].cub[selection[1][0]].read_file(true, false, false);
							wavy[selection[0][0]].cub[selection[1][0]].write_dgrid(wavy[selection[0][0]], path, opt.debug);
							break;
						}
						case 'X':
						case 'x': {
							if (wavy[activewave].cub.size() >= 1) {
								std::cout << "Which cube do you want to save?" << endl;
								selection.resize(2);
								selection[0].resize(1);
								selection[1].resize(1);
								select_cubes(selection, wavy, 1);
							}
							string path = wavy[selection[0][0]].cub[selection[1][0]].path;
							path = path.substr(0, path.find(".cub"));
							path += ".xdgrid";
							if (!wavy[selection[0][0]].cub[selection[1][0]].get_loaded())
								wavy[selection[0][0]].cub[selection[1][0]].read_file(true, false, false);
							wavy[selection[0][0]].cub[selection[1][0]].write_xdgraph(wavy[selection[0][0]], path, opt.debug);
							break;
						}
						default:
							std::cout << "Sorry, i did not get that! Try again!" << endl;
							continue;
						}
						break;
					}
				}
				else {
					if ((wavy[activewave].get_origin() == 2 || wavy[activewave].get_origin() == 4) && w_or_c) {
						std::cout << "Which format do you want to save the wavefunction in?" << endl
							<< "W) WFN format" << endl
							<< "F) Fchk format" << endl;
						std::cin >> sel;
						switch (sel) {
						case 'W':
						case 'w': {
							endings.push_back(".wfn");
							endings.push_back(".ffn");
							string path;
							if (!expert) save_file_dialog(path, opt.debug, endings);
							else {
								std::cout << "Enter filename: ";
								std::cin >> path;
								while (exists(path)) {
									std::cout << path << " exists, do you want to overwrite it? ";
									if (!yesno()) {
										std::cout << "Then try again: ";
										std::cin >> path;
									}
								}
							}
							bool all = false;
							if (path.find(".wfn") != -1) {
								std::cout << "Do you want to write all MOs?" << endl;
								all = yesno();
							}
							if (!writewfn(wavy[activewave], path, opt.debug, !all)) {
								Enter();
								cls();
							}
							else {
								if (opt.debug) Enter();
								cls();
								std::cout << "Wrote Wavefunction!\n";
							}
							break;
						}
						case 'F':
						case 'f': {
							//opt.debug=true;
							endings.push_back(".fchk");
							endings.push_back(".Fchk");
							endings.push_back(".FChk");
							string outputname = wavy[activewave].get_path();
							if (opt.debug)std::cout << "Loaded path..." << endl;
							size_t where;
							if (wavy[activewave].get_origin() == 2) where = outputname.find("wfn");
							else if (wavy[activewave].get_origin() == 1) where = outputname.find("out");
							else if (wavy[activewave].get_origin() == 4) where = outputname.find("ffn");
							if (opt.debug)std::cout << "Found the file ending in: " << outputname << " at position: " << where << "; npos= " << string::npos << " Origin is: " << wavy[activewave].get_origin() << endl;
							if (where >= outputname.length()) {
								if (!expert) save_file_dialog(outputname, opt.debug, endings);
								else {
									std::cout << "Enter filename: ";
									std::cin >> outputname;
									while (exists(outputname)) {
										std::cout << outputname << " exists, do you want to overwrite it? ";
										if (!yesno()) {
											std::cout << "Then try again: ";
											std::cin >> outputname;
										}
									}
								}
							}
							else outputname.erase(where, 3);
							string basis_temp = wavy[activewave].get_basis_set_name();
							if (basis_temp.length() < 3) {
								int tries = 0;
								bool end = false;
								string temp;
								while (!end && tries != 3) {
									std::cout << "Please give the name of the basis set you want to use: ";
									std::cin >> temp;
									string basis_set_file(basis_set_path);
									basis_set_file.append(temp);
									if (opt.debug)std::cout << "looking for: " << basis_set_file << endl;
									if (exists(basis_set_file)) end = true;
									else tries++;
								}
								if (tries == 3) {
									std::cout << "Sorry, this takes too long... please make sure you know what you want and try again!" << endl;
									Enter();
									break;
								}
								wavy[activewave].change_basis_set_name(temp);
							}
							bool read = false;
							if (expert) {
								std::cout << "What is the charge of your molecule?" << endl;
								int temp = 0;
								std::cin >> temp;
								wavy[activewave].assign_charge(temp);
							}
							else wavy[activewave].assign_charge(wavy[activewave].calculate_charge());
							if (wavy[activewave].get_multi() == 0) wavy[activewave].guess_multiplicity(std::cout);
							free_fchk(std::cout, outputname, opt.basis_set_path, wavy[activewave], opt.debug);
							break;
						}
						default: {
							std::cout << "Sorry, i didn't get that, could you try it again?\n";
							Enter();
							cls();
						}
						}
					}
					else if (wavy[activewave].get_origin() == 3 || !w_or_c) {
						string path;
						if (wavy[activewave].cub.size() <= 0) {
							std::cout << "No cubes loaded!" << endl;
							break;
						}
						int nr1 = 0;
						if (wavy[activewave].cub.size() >= 1) {
							std::cout << "Which cube do you want to save?" << endl;
							selection.resize(2);
							selection[0].resize(1);
							selection[1].resize(1);
							if (!expert) select_cubes(selection, wavy, 1);
							else {
								for (int i = 0; i < wavy[selection[0][0]].cub.size(); i++) {
									std::cout << " " << i << ") " << wavy[activewave].cub[i].path;
									if (!exists(wavy[activewave].cub[i].path))std::cout << " (MEM ONLY)";
									std::cout << endl;
								}
								std::cout << "Cube 1: "; cin >> nr1;
								while (nr1 >= wavy[activewave].cub.size() || nr1 < 0) {
									std::cout << "Invalid choice, select again: ";
									std::cin >> nr1;
								}
								selection[0][0] = activewave;
								selection[1][0] = nr1;
							}
						}

						endings.push_back(".cube");
						endings.push_back(".cub");
						if (!expert)
							save_file_dialog(path, opt.debug, endings);
						else {
							std::cout << "Give filepath please: ";
							std::cin >> path;
						}
						wavy[selection[0][0]].cub[selection[1][0]].write_file(wavy[selection[0][0]], path, opt.debug);
					}
				}
				endings.resize(0);
				break;
			}
			case 'O':
			case 'o':{
				if(wavy.size()<1){
					std::cout << "First you need to read a wavefunction!" << endl;
					break;
				}
				if(wavy[activewave].get_origin ()==2||wavy[activewave].get_origin()==4){
					std::cout << "Sorting wavefunction!" << endl;
					wavy[activewave].sort_wfn (wavy[activewave].check_order (opt.debug),opt.debug);
					if(opt.debug) Enter();
					cls();
				}
				else{
					std::cout << "I can only sort .wfn/ffn files!" << endl;
					Enter();
					std::cout << endl;
				}
				break;
			}
			case 'B':
			case 'b':{
				if(wavy.size()<1){
					std::cout << "First you need to read a wavefunction!" << endl;
					break;
				}
				if(!read_basis_set_vanilla (opt.basis_set_path,wavy[activewave],opt.debug))std::cout << "Problem during reading of the basis set!" << endl;
				if(opt.debug) Enter();
				cls();
				break;
			case 'a':
			case 'A':
				if(wavy.size()<1) {
					cls();
					break;
				}
				vector < vector < unsigned int > > selection;
				selection.resize(1);
				select_cubes(selection,wavy,1,true);
				activewave=selection[0][0];
				cls();
				break;
			}
			case 'e':
			case 'E':
				if(expert) expert=false;
				else expert=true;
				break;
			case 'q':
			case 'Q':
				if(!unsaved_files(wavy))
					end=true;
				else{
					std::cout << "There are unsaved files! Are you sure you want to exit? ";
					if(yesno())
						end=true;
				}
				break;
			case 'd':
			case 'D':
				opt.debug=true;
				cls();
				break;
			case 'u':
			case 'U':
				if(check_bohr (wavy[activewave],opt.debug)){
					cls();
					std::cout << "Appears to be in bohr!" << endl;
				}
				else {
					cls();
					std::cout << "Appears to be in Angström!" << endl;
				}
				break;
			default:
				if(opt.debug){
					std::cout << "This command is unknown!" << endl;
					Enter();
				}
				cls();
				std::cout << "Sorry, i didn't get that, could you try it again?\n";
				break;
		}			
	}
	cls();
	return 0;
}

