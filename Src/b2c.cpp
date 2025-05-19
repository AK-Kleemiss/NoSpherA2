#include "pch.h"
#include "wfn_class.h"
#include "convenience.h"

using namespace std;

double vec_length(vector <double> &vector){        //    Calculates the length of a given Vector (obviously...)
    return sqrt(abs(vector[0]*vector[0]
            +vector[1]*vector[1]
            +vector[2]*vector[2]));
}

struct cubepoint{
    int x;
    int y;
    int z;
    double value;
};

bool b2c(const cube* cub, const vector<atom> &atoms, bool debug, bool bcp)
{
    iMatrix3 CP(cub->get_size(0), cub->get_size(1), cub->get_size(2));
    ivec2 Liste(3);
    double GradMax, xlength, ylength, zlength;
    xlength = sqrt(cub->get_vector(0, 0)*cub->get_vector(0, 0) + cub->get_vector(0, 1)*cub->get_vector(0, 1) + cub->get_vector(0, 2)*cub->get_vector(0, 2));
    ylength = sqrt(cub->get_vector(1, 0)*cub->get_vector(1, 0) + cub->get_vector(1, 1)*cub->get_vector(1, 1) + cub->get_vector(1, 2)*cub->get_vector(1, 2));
    zlength = sqrt(cub->get_vector(2, 0)*cub->get_vector(2, 0) + cub->get_vector(2, 1)*cub->get_vector(2, 1) + cub->get_vector(2, 2)*cub->get_vector(2, 2));
    if(debug)std::cout << "calculated lengths!" << endl;
    string s, s1, s2, s3;
    double xmin = cub->get_origin(0);//, xstep;
    double ymin = cub->get_origin(1);//, ystep;
    double zmin = cub->get_origin(2);//, zstep;

    if(debug)std::cout << "Resized CP and initialized it; GOING TO MAKE LISTE NOW" << endl;
    int iCP=0, ListeMax;
    cubepoint Max={0,0,0,0.0};
    vec Maxima;
    vec3 distances;
    svec Labels;
    ivec nrs;
    vector<vector<bvec>> border;
    ivec2 neighbours;
    vector<vector<cubepoint> > BCPs;
    distances.resize(3);
    for(int i=0; i<3; i++){
        distances[i].resize(3);
        for(int j=0; j<3; j++){
            distances[i][j].resize(3);
            for(int k=0; k<3; k++) distances[i][j][k]=sqrt(
                     pow((1-i)*xlength,2)
                    +pow((1-j)*ylength,2)
                    +pow((1-k)*zlength,2));
        }
    }
    if(debug)std::cout << "calculated distances!" << endl;
    for (int x=0; x<cub->get_size(0); x++)
        for (int y=0; y<cub->get_size(1); y++)
            for (int z=0; z<cub->get_size(2); z++){
                if (CP(x,y,z) == 0){
                    ListeMax = 0;
                    Max.x = x; Max.y = y; Max.z = z;
                    int xs, ys, zs;
                    do{
                        ListeMax++;
                        xs=Max.x, ys=Max.y, zs=Max.z;
                        if(xs<cub->get_size(0) && ys<cub->get_size(1) && zs<cub->get_size(2)) 
                            Max.value=(cub->get_value(xs,ys,zs));
                        else 
                            Max.value=0.0;
                        if(Liste[0].size()<=ListeMax) 
                            Liste[0].push_back(Max.x);  //only append List when needed check if neccesarry
                        else 
                            Liste[0][ListeMax-1]=Max.x;
                        if(Liste[1].size()<=ListeMax) 
                            Liste[1].push_back(Max.y);
                        else 
                            Liste[1][ListeMax-1]=Max.y;
                        if(Liste[2].size()<=ListeMax) 
                            Liste[2].push_back(Max.z);
                        else 
                            Liste[2][ListeMax-1]=Max.z;
                        GradMax=0;
                        for (int ix=xs-1; ix<xs+2; ix++)
                            for (int iy=ys-1; iy<ys+2; iy++)
                                for (int iz=zs-1; iz<zs+2; iz++){
                                    if (ix==xs && iy==ys && iz==zs) continue;
                                    if (ix<0 || iy<0 || iz<0) continue;
                                    if (ix>=cub->get_size(0) || iy>=cub->get_size(1) || iz>=cub->get_size(2)) continue;
                                    /*double dist=sqrt(
                                             ((ix-xs)*xlength)*((ix-xs)*xlength)
                                            +((iy-ys)*ylength)*((iy-ys)*ylength)
                                            +((iz-zs)*zlength)*((iz-zs)*zlength));*/
                                    //    Tests if this voxel has maximum value compared to neighboring voxels
                                    if (((cub->get_value(ix, iy, iz) - (cub->get_value(xs, ys, zs))) / distances[1 + ix - xs][1 + iy - ys][1 + iz - zs]) > GradMax) {
                                        Max.x = ix; Max.y = iy; Max.z = iz;
                                        Max.value = cub->get_value(ix, iy, iz);
                                        GradMax = ((cub->get_value(ix, iy, iz)) - (cub->get_value(xs, ys, zs))) / distances[1 + ix - xs][1 + iy - ys][1 + iz - zs];
                                    }
                                }
                    } while (!(GradMax == 0 || CP(Max.x, Max.y, Max.z) > 0));
                    if (CP(Max.x, Max.y, Max.z) > 0) for (int i = 0; i < ListeMax; i++) CP(Liste[0][i], Liste[1][i], Liste[2][i]) = CP(Max.x, Max.y, Max.z);
                    else{
                        double pos[3]={Max.x*xlength+cub->get_origin(0),Max.y*ylength+cub->get_origin(1),Max.z*zlength+cub->get_origin(2)};
                        if (debug) {
                            std::cout << "DBUG: Position of CP: ";
                            for (int i=0; i<3; i++)std::cout << pos[i] << " ";
                            std::cout << endl;
                        }
                        double min_dist, temp;
                        unsigned int atom = 0;
                        Maxima.push_back(Max.value);
                        iCP++;
                        min_dist = 10000;
                        for (int i = 0; i < ListeMax; i++) CP(Liste[0][i], Liste[1][i], Liste[2][i]) = iCP;
                        for (int i=0; i<atoms.size(); i++){
                            temp=sqrt(pow(pos[0]-atoms[i].get_coordinate(0), 2) + pow(pos[1] - atoms[i].get_coordinate(1), 2) + pow(pos[2] - atoms[i].get_coordinate(2), 2));
                            //if(debug)std::cout << "DBUG: Distance to atom " << i  << " with label " << atoms[i].get_label() << " is " << temp;
                            if(min_dist>temp){
                                min_dist=temp;
                                atom=i;
                                //if(debug)std::cout << " *";
                            }
                            //if(debug)std::cout << endl;
                        }
                        Labels.push_back(atoms[atom].get_label());
                        nrs.push_back(atom);
                    }
                    for(int i=0; i<3; i++) 
                        Liste[i].resize(0); //Revert Liste and free memory
                }
            }
    if(bcp){
        neighbours.resize(iCP);
        BCPs.resize(iCP);
        border.resize(cub->get_size(0));
        for(int i=0; i<cub->get_size(0); i++){
            border[i].resize(cub->get_size(1));
            for(int j=0; j<cub->get_size(1); j++){
                border[i][j].resize(cub->get_size(2));
                for(int k=0; k<cub->get_size(2);k++) 
                    border[i][j][k]=false;
            }
        }
        for(int x=0; x<cub->get_size(0); x++)
            for(int y=0; y<cub->get_size(1); y++)
                for(int z=0; z<cub->get_size(2); z++)
                    for (int ix=x-1; ix<x+2; ix++)
                        for (int iy=y-1; iy<y+2; iy++)
                            for (int iz=z-1; iz<z+2; iz++){
                                if (ix==x && iy==y && iz==z) continue;
                                if (ix<=0 || iy<=0 || iz<=0) continue;
                                if (ix>=cub->get_size(0) || iy>=cub->get_size(1) || iz>=cub->get_size(2)) continue;
                                //    Tests if this voxel has neighboring voxels, that are from a different basin and assigns list of neighbors for each basin
                                if (CP(ix, iy, iz) != CP(x, y, z)){
                                    border[x][y][z]=true;
                                    int found=false;
                                    for (int i = 0; i < neighbours[CP(x, y, z)].size(); i++) 
                                        if (neighbours[CP(x, y, z)][i] == CP(ix, iy, iz)) 
                                            found = true;
                                    if (!found) 
                                        neighbours[CP(x, y, z)].push_back(CP(ix, iy, iz));
                                }
                            }
        //sanity check, all basins must be neighboring each other pairwise!
        for(int b=0; b<iCP; b++)
            for(int n=0; n<neighbours[b].size(); n++)
                if(find(neighbours[n].begin(), neighbours[n].end(),b)==neighbours[n].end()){
                    std::cout << "ERROR: Basins should be neighbours in a pairwise way! Basin " << b << " has neighbour " << n << ", but it does not appear to be the case the other way around!" << endl;
                    return false;
                }
        //possibly better if we assigned max and min values for coordinates of basins to reduce amount of calculations, lets see...
        for(int b=0; b<iCP; b++){
            Max={0,0,0,0.0};
            for(int n=0; n<neighbours[b].size(); n++){
                for(int x=0; x<cub->get_size(0); x++)
                    for(int y=0; y<cub->get_size(1); y++)
                        for(int z=0; z<cub->get_size(2); z++){
                            //check that we are on A) A border B) in the basin we want C) the neighboring basin we look at is at this border D) the value is bigger than the previously found maximum
                            if(!border[x][y][z]) 
                                continue;
                            if(CP(x, y, z) !=b) 
                                continue;
                            bool found=false;
                            for (int ix=x-1; ix<x+2; ix++)
                                for (int iy=y-1; iy<y+2; iy++)
                                    for (int iz=z-1; iz<z+2; iz++){
                                        if (ix==x && iy==y && iz==z) continue;
                                        if (ix<=0 || iy<=0 || iz<=0) continue;
                                        if (ix>=cub->get_size(0) || iy>=cub->get_size(1) || iz>=cub->get_size(2)) continue;
                                        if (CP(ix, iy, iz) !=n) continue;
                                        else found=true;
                                    }
                            if(found&&cub->get_value(x,y,z)>Max.value) 
                                Max={x,y,z,cub->get_value(x,y,z)};
                        }
                BCPs[b].push_back(Max);
            }
        }
        for(int b=0; b<iCP; b++)
            for(int n=0; n<neighbours[b].size(); n++){
                int back_reference=0;
                for (int i = 0; i < neighbours[neighbours[b][n]].size(); i++)
                    if (neighbours[neighbours[b][n]][i] == b)
                        back_reference = i;
                if(BCPs[b][n].x-BCPs[neighbours[b][n]][back_reference].x>1 || BCPs[b][n].y-BCPs[neighbours[b][n]][back_reference].y > 2 || BCPs[b][n].z-BCPs[neighbours[b][n]][back_reference].z > 2){
                    std::cout << "The BCPs on both sides of the ZFS are too far apart..." << endl;
                }
                else{
                    std::cout << "BCP: " << Labels[b] << "-" << Labels[neighbours[b][n]] << " ED: " << (BCPs[b][n].value + BCPs[neighbours[b][n]][back_reference].value) / 2 << " Position: "
                        << (BCPs[b][n].x + BCPs[neighbours[b][n]][back_reference].x) / 2 * cub->get_vector(0, 0)
                        + (BCPs[b][n].y + BCPs[neighbours[b][n]][back_reference].y) / 2 * cub->get_vector(1, 0)
                        + (BCPs[b][n].z + BCPs[neighbours[b][n]][back_reference].z) / 2 * cub->get_vector(2, 0)
                        + xmin << " "
                        << (BCPs[b][n].x + BCPs[neighbours[b][n]][back_reference].x) / 2 * cub->get_vector(0, 1)
                        + (BCPs[b][n].y + BCPs[neighbours[b][n]][back_reference].y) / 2 * cub->get_vector(1, 1)
                        + (BCPs[b][n].z + BCPs[neighbours[b][n]][back_reference].z) / 2 * cub->get_vector(2, 1)
                        + ymin
                        << (BCPs[b][n].x + BCPs[neighbours[b][n]][back_reference].x) / 2 * cub->get_vector(0, 2)
                        + (BCPs[b][n].y + BCPs[neighbours[b][n]][back_reference].y) / 2 * cub->get_vector(1, 2)
                        + (BCPs[b][n].z + BCPs[neighbours[b][n]][back_reference].z) / 2 * cub->get_vector(2, 2)
                        + zmin
                        << endl;
                }
            }
        if(debug)std::cout << "done with BCPs" << endl;
    }
    if(debug)std::cout << "done with liste, writing basins now!" << endl;
    ivec basins;
    vec EDS(iCP);
    vec VOL(iCP);
    double dv = cub->get_dv();
    if(debug)std::cout << "dv: " << dv  << " iCP: " << iCP << endl;
    for (int x=0; x<cub->get_size(0); x++)
        for (int y=0; y<cub->get_size(1); y++)
            for (int z=0; z<cub->get_size(2); z++)
                for (int a = 0; a < iCP; a++) if (a + 1 == CP(x, y, z)) {
                    EDS[a] += cub->get_value(x, y, z) * dv;            //    Electrons in this voxel
                    VOL[a] += dv;                                    //    Size of the voxel
                }
    if(debug) for(int a=0; a<iCP; a++)std::cout << a << " eds: " << EDS[a] << " vol: " << VOL[a] << endl;
    unsigned int in=9999999;
    std::cout << "I found " << iCP << " Basins." << endl;
    for(int a=0; a<iCP; a++)std::cout << "Integrated value in Basin " << toString<int>(a+1) << ": " << scientific << setw(14) << setprecision(7) << EDS[a]
                                  << " Volume: " << VOL[a] << " Maximum: " << Maxima[a] << " possibly atom: " << Labels[a] << "_" << nrs[a]<< endl;
    if(debug) {
        std::cout << "DEBUG: Labels_size: " << Labels.size() << endl;
        for(int i=0; i< Labels.size(); i++)std::cout << Labels[i] << " ";
        std::cout << endl;
    }
    std::cout << "Which of these do you want to include into one cube file?" << endl
         << "Put 0 for end of input. If none are selected (first number is zero), all basins will be written into separate files" << endl;
    while (in!=0){
        cin >> in;
        if(in!=0) basins.push_back(in);
    }
    string temp;
    string replace("E");
    string path_temp=cub->get_path().replace_extension(".b2c_log").generic_string();
    path_temp.erase(path_temp.find(".cub"),5);
    ofstream logfile(path_temp.c_str(),ios::out);
    logfile << "Number of Basins: " << iCP << endl;
    double Integral=0.0;
    temp=cub->get_path().generic_string();
    temp.erase(temp.find(".cub"),5);
    if(basins.size()>0){
        temp+="_"+toString<int>((int) basins.size())+"_basins";
        //for(int i=0; i<basins.size(); i++){
        //    temp+='_'+Labels[basins[i]]+'_'+toString<unsigned int>(nrs[basins[i]-1]);
        //   std::cout << temp << endl;
        //}
        temp+=".cube";
        ofstream outfile(temp.c_str(), ios::out);
        outfile << "First comment line?" << endl << "second comment line" << endl << "   " << atoms.size();
        outfile << fixed;
        outfile << setprecision(6);
        outfile << setw(12);
        outfile << xmin << " " << ymin << " " << zmin << endl; // output comments and origin
        for(int i=0; i<3; i++){
            outfile << setprecision(0);
            outfile << setw(6);
            outfile << cub->get_size(i) << " ";
            outfile << setw(12);
            outfile << setprecision (6);
            for(int j=0; j<3; j++) outfile << cub->get_vector(i,j) << " "; //    Output principal axis and size of cube
            outfile << endl;
            outfile.flush();
        }
        for (int i=0; i<atoms.size(); i++){
            outfile << "  " << atoms[i].get_charge();
            outfile << setw(12) << setprecision(6);
            outfile << atoms[i].get_charge() << " ";
            outfile << atoms[i].get_coordinate(0) << " " << atoms[i].get_coordinate(1) << " " << atoms[i].get_coordinate(2) << " ";    //    write atoms
            outfile << endl;
        }
        outfile.flush();
        stringstream stream;
        for (int x=0; x<cub->get_size(0); x++){
            for (int y=0; y<cub->get_size(1); y++){
                   unsigned int r=0;
                   for (int z=0; z<cub->get_size(2); z++){
                       bool include=false;
                    for (int a = 0; a < basins.size(); a++) 
                        if (basins[a] == CP(x, y, z)) 
                            include = true;
                    stream << uppercase << scientific << setw(14) << setprecision(7) << cub->get_value(x, y, z) * include;
                    outfile << stream.str();
                    stream.str("");
                    r++;
                    if(r%6==0) outfile << endl;
                   }
                   if(r%6!=0) outfile << endl;
            }
        outfile.flush();
        }
        for(int a=0; a<basins.size(); a++){
            logfile << "Integrated value in Basin " << toString<int>(basins[a]) << ": " << scientific << setw(14) << setprecision(7) << EDS[basins[a]-1]
                << "Volume: " << VOL[basins[a]-1] << " Atom: " << Labels[basins[a]-1] << "_" << nrs[basins[a]-1] << endl;
            Integral+=EDS[basins[a]-1];
        }
        outfile.close();
    }
    else{
        for(int f=0; f<iCP; f++){
            temp=cub->get_path().generic_string() + '_' + Labels[f] + '_' + toString<unsigned int>(nrs[f]) + '_' + toString<int>(f) + ".cube";
            temp.erase(temp.find(".cub"),5);
            ofstream outfile(temp.c_str(), ios::out);
            outfile << s1 << "\n" << s2 << "\n" << setw(5) << atoms.size();
            outfile << fixed << setprecision(6) << setw(12) << xmin << " " << ymin << " " << zmin << "\n"; // output comments and origin
            for(int i=0; i<3; i++)
                outfile << setprecision(0) << setw(5) << cub->get_size(i) << setw(12) << setprecision (6) << cub->get_vector(i,0) << setw(12) << setprecision(6) << cub->get_vector(i, 1) << setw(12) << setprecision(6) << cub->get_vector(i, 2) << "\n"; //    Output principal axis and size of cube
            for (int i=0; i<atoms.size(); i++){
                outfile << setw(5) << atoms[i].get_charge() << setw(5) << atoms[i].get_charge() << ".000000";
                outfile << fixed << setw(12) << setprecision(6) << atoms[i].get_coordinate(0);
                outfile << fixed << setw(12) << setprecision(6) << atoms[i].get_coordinate(1);
                outfile << fixed << setw(12) << setprecision(6) << atoms[i].get_coordinate(2);
                outfile << "\n";
            }
            for (int x=0; x<cub->get_size(0); x++){
                for (int y=0; y<cub->get_size(1); y++){
                    unsigned int r=0;
                    for (int z=0; z<cub->get_size(2); z++){
                        // Write values into cube
                        outfile << uppercase << scientific << setw(13) << setprecision(5) << cub->get_value(x, y, z) * (CP(x, y, z) == (f + 1));
                        r++;
                        if(r%6==0) outfile << "\n";
                    }
                    if(r%6!=0) outfile << "\n";
                }
            }
            outfile.flush();
            logfile << "Integrated value in Basin " << toString<int>(f+1) << ": " << scientific << setw(14) << setprecision(7) << EDS[f] << "Volume: " << VOL[f]
                    << " Atom: " << Labels[f] << "_" << nrs[f] << endl;
            outfile.close();
            Integral+=EDS[f];
        }
    }
    logfile << "Integral over all basins: " << Integral << endl;
    if (bcp)
        for (int b=0; b<iCP; b++)
            for (int n=0; n<neighbours[b].size(); n++){
                int back_reference=0;
                for (int i=0; i<neighbours[neighbours[b][n]].size(); i++) if(neighbours[neighbours[b][n]][i]==b) back_reference=i;
                logfile << "BCP: " << Labels[b] << "-" << Labels[neighbours[b][n]] << " ED: " << (BCPs[b][n].value+BCPs[neighbours[b][n]][back_reference].value)/2  << " Position: "
                        << (BCPs[b][n].x+BCPs[neighbours[b][n]][back_reference].x)/2*cub->get_vector(0,0)
                         + (BCPs[b][n].y+BCPs[neighbours[b][n]][back_reference].y)/2*cub->get_vector(1,0)
                         + (BCPs[b][n].z+BCPs[neighbours[b][n]][back_reference].z)/2*cub->get_vector(2,0)
                         + xmin << " "
                        << (BCPs[b][n].x+BCPs[neighbours[b][n]][back_reference].x)/2*cub->get_vector(0,1)
                         + (BCPs[b][n].y+BCPs[neighbours[b][n]][back_reference].y)/2*cub->get_vector(1,1)
                         + (BCPs[b][n].z+BCPs[neighbours[b][n]][back_reference].z)/2*cub->get_vector(2,1)
                         + ymin
                        << (BCPs[b][n].x+BCPs[neighbours[b][n]][back_reference].x)/2*cub->get_vector(0,2)
                         + (BCPs[b][n].y+BCPs[neighbours[b][n]][back_reference].y)/2*cub->get_vector(1,2)
                         + (BCPs[b][n].z+BCPs[neighbours[b][n]][back_reference].z)/2*cub->get_vector(2,2)
                         + zmin
                        << endl;
            }
    logfile.flush();
    logfile.close();
    return true;
};
