/**********************
* main.cpp
* Jeongho Son
***********************/

#define TRUE 1
#define FALSE 0

#include <vcg/complex/complex.h>
#include <vcg/simplex/face/topology.h>
#include <vcg/simplex/face/jumping_pos.h>

#include <QString>
#include <QDir>

#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>

#include <iostream>
#include <string>
#include <cstdlib>

class MyVertex;
class MyEdge;
class MyFace;

struct MyUsedTypes : public vcg::UsedTypes<	vcg::Use<MyVertex>   ::AsVertexType,
                                        vcg::Use<MyEdge>     ::AsEdgeType,
                                        vcg::Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public vcg::Vertex<MyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::VFAdj,vcg::vertex::BitFlags, vcg::vertex::Mark  >{};
class MyFace    : public vcg::Face< MyUsedTypes, vcg::face::VFAdj, vcg::face::FFAdj,  vcg::face::VertexRef, vcg::face::BitFlags > {};
class MyEdge    : public vcg::Edge<MyUsedTypes>{};
class MyMesh    : public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyEdge>,std::vector<MyFace>  > {};

typedef typename MyMesh::CoordType CoordType;

static int print_verbose = 0;
static int p_level = 4;
static int prc_embed = -1;
static int prc_extract = -1;
static char *ori_mesh_name = NULL;
static char *stego_mesh_name = NULL;
static char *recov_mesh_name = NULL;
static char *wm_file_name = NULL;

void OpenMesh( MyMesh &m, char *file_name );

void PrintSNR( MyMesh &m, MyMesh &n )
{
    MyMesh::VertexIterator mvi;
    MyMesh::VertexIterator nvi;
    
    double mx = 0, my = 0, mz = 0;
    double nx = 0, ny = 0, nz = 0;
    double mc = 0, nm = 0;
    for(mvi=m.vert.begin();mvi!=m.vert.end();++mvi){
        mx += (*mvi).cP()[0]; my += (*mvi).cP()[1]; mz += (*mvi).cP()[2];
    }
    mx /= m.VN(); my /= m.VN(); mz /= m.VN();
    for(nvi=n.vert.begin();nvi!=n.vert.end();++nvi){
        nx += (*nvi).cP()[0]; ny += (*nvi).cP()[1]; nz += (*nvi).cP()[2];
    }
    nx /= n.VN(); ny /= n.VN(); nz /= n.VN();
    for(mvi=m.vert.begin();mvi!=m.vert.end();++mvi){
        mc += std::pow((*mvi).cP()[0] - mx, 2) + std::pow((*mvi).cP()[1] - my, 2) + std::pow((*mvi).cP()[2] - mz, 2);
    }
    mc /= m.VN();
    nvi = n.vert.begin();
    for(mvi=m.vert.begin();mvi!=m.vert.end();++mvi){
        nm += std::pow((*mvi).cP()[0] - (*nvi).cP()[0], 2) + std::pow((*mvi).cP()[1] - (*nvi).cP()[1], 2) + std::pow((*mvi).cP()[2] - (*nvi).cP()[2], 2);
        nvi++;
    }
    nm /= m.VN();
    std::cout << "SNR: " << mc/nm << std::endl;
}

double CalculateDist( CoordType x, CoordType y){
    double dist;
    dist = sqrt(pow((x[0]-y[0]),2) + pow((x[1]-y[1]),2)+pow((x[2]-y[2]),2));
    return dist;
}


void _StartProcess( MyMesh &m )
{
    // embed wm on m
    // finding unit EU_i based on the connectivity
    vcg::tri::UpdateTopology<MyMesh>::FaceFace(m);
    vcg::tri::UpdateTopology<MyMesh>::VertexFace(m);
    MyMesh::VertexIterator vi;
    int checkflag;
    int i=0;
    int deg_vi;
    //index marking
    
    for(vi=m.vert.begin();vi!=m.vert.end();++vi){
        i++; (*vi).IMark() = i;
        //std::cout << (*vi).cIMark() << std::endl;
        //std::cout << (*vi).cFlags() << std::endl;        
    }
    vi = m.vert.begin(); 
    
    while(vi!=m.vert.end()){
        //if(!vcg::tri::io::ExporterPLY<MyMesh>::Save(m, "temp3.ply")) {std::cout << (*vi).cIMark() << std::endl;}
        //std::cout << (*vi).cIMark() << std::endl;
        if((*vi).cFlags()) {vi++; continue;}
        //neighborhood check
        MyVertex* cvi;
        MyVertex* fNvi;
        MyVertex* Nvi;
        cvi = &(*vi);
        vcg::face::JumpingPos<MyFace> pos((*vi).VFp(),cvi);
        fNvi = pos.VFlip(); Nvi = fNvi;
        checkflag = 0;
        do {
            checkflag += (*Nvi).cFlags();
            pos.FlipF();
            pos.FlipE();
            Nvi = pos.VFlip();
        } while(Nvi != fNvi);
        if(checkflag) {vi++; continue;}
           
        // do embedding
        
        vcg::face::JumpingPos<MyFace> pos1((*vi).VFp(),cvi);
        fNvi = pos1.VFlip();
        Nvi = fNvi;
        CoordType pred(0,0,0);
        CoordType emb_vi(0,0,0);
        (*vi).Flags()=1;
        deg_vi = 0;
        do {
            deg_vi++;
            pred = pred + (*Nvi).cP();
            (*Nvi).Flags() = 1;
            pos1.FlipF();
            pos1.FlipE();
            Nvi = pos1.VFlip();
        } while(Nvi != fNvi);
        pred = pred/deg_vi;
        
        //std::cout << pred[1] << "  " << (*vi).P()[1] << std::endl;
        double di = CalculateDist(pred, (*vi).cP())/2;
        //std::cout << di << std::endl;
        double di_t = di + ((int)(di*std::pow(10,p_level)))/std::pow(10,p_level) + 1/(std::pow(10,p_level));
        emb_vi = pred + ((*vi).cP() - pred) * (di_t/(2*di));

        //std::cout << "change " << (*vi).cIMark() << std::endl;
        //std::cout << (*vi).cP()[0] << (*vi).cP()[1] << (*vi).cP()[2] << std::endl;
        //std::cout << emb_vi[0] << emb_vi[1] << emb_vi[2] << std::endl;

        (*vi).P() = emb_vi;
        
        vi++;
        
    }
    for(vi=m.vert.begin();vi!=m.vert.end();++vi){
        (*vi).Flags() = 0;
        //i++; (*vi).IMark() = i;
        //std::cout << (*vi).cIMark() << std::endl;
        //std::cout << (*vi).cFlags() << std::endl; 
        //std::cout << (*vi).cP()[0] << (*vi).cP()[1] << (*vi).cP()[2] << std::endl;       
    }
    std::cout << m.VN() << std::endl;
    //std::cout << i << std::endl;
    vcg::tri::io::Exporter<MyMesh>::Save(m, "temp.ply");

}

void StartProcess()
{
    std::cout << "start processing.. " << std::endl;
    if(prc_embed){
        // start embedding on (ori) and write to (stego) using (wm, p_level)
        
        // io check
        if( !ori_mesh_name || !stego_mesh_name || !wm_file_name ) {
            std::cerr << "need -ori ori_mesh_name -stego stego_mesh_name -wm wm_file_name" << std::endl;
            return;
        }
        MyMesh m;
        OpenMesh(m, ori_mesh_name);
        
        vcg::tri::UpdateTopology<MyMesh>::FaceFace(m);
        vcg::tri::UpdateTopology<MyMesh>::VertexFace(m);
        MyMesh::VertexIterator vi;
        int checkflag;
        int i=0;
        int deg_vi;

        //index marking
        for(vi=m.vert.begin();vi!=m.vert.end();++vi){
            i++; (*vi).IMark() = i;
        }
        vi = m.vert.begin(); 
        
        while(vi!=m.vert.end()){

            if((*vi).cFlags()) {vi++; continue;}
            //neighborhood check
            MyVertex* cvi;
            MyVertex* fNvi;
            MyVertex* Nvi;
            cvi = &(*vi);
            vcg::face::JumpingPos<MyFace> pos((*vi).VFp(),cvi);
            fNvi = pos.VFlip(); Nvi = fNvi;
            checkflag = 0;
            do {
                checkflag += (*Nvi).cFlags();
                pos.FlipF();
                pos.FlipE();
                Nvi = pos.VFlip();
            } while(Nvi != fNvi);
            if(checkflag) {vi++; continue;}
               
            // do embedding
            
            vcg::face::JumpingPos<MyFace> pos1((*vi).VFp(),cvi);
            fNvi = pos1.VFlip();
            Nvi = fNvi;
            CoordType pred(0,0,0);
            CoordType emb_vi(0,0,0);
            (*vi).Flags()=1;
            deg_vi = 0;
            do {
                deg_vi++;
                pred = pred + (*Nvi).cP();
                (*Nvi).Flags() = 1;
                pos1.FlipF();
                pos1.FlipE();
                Nvi = pos1.VFlip();
            } while(Nvi != fNvi);
            pred = pred/deg_vi;
            
            //std::cout << pred[1] << "  " << (*vi).P()[1] << std::endl;
            double di = CalculateDist(pred, (*vi).cP())/2;
            //std::cout << di << std::endl;
            double di_t = di + ((int)(di*std::pow(10,p_level)))/std::pow(10,p_level) + 1/(std::pow(10,p_level));
            emb_vi = pred + ((*vi).cP() - pred) * (di_t/(1*di));
            (*vi).P() = emb_vi;
            
            vi++;
            
        }
        for(vi=m.vert.begin();vi!=m.vert.end();++vi){
            (*vi).Flags() = 0;       
        }
        std::cout << m.VN() << std::endl;
        //std::cout << i << std::endl;
        vcg::tri::io::Exporter<MyMesh>::Save(m, stego_mesh_name);
        
    }
    else if(prc_extract){
        // start extracting on (stego) and write to (recov) using (p_level)
        // if needed, start analyzing (wm_ext,recov) with (wm,ori) to evaluate the performance
        if( !recov_mesh_name || !stego_mesh_name || !wm_file_name ) {
            std::cerr << "need -recov recov_mesh_name -stego stego_mesh_name -wm wm_file_name" << std::endl;
            return;
        }
        MyMesh m;
        OpenMesh(m, stego_mesh_name);
        
        vcg::tri::UpdateTopology<MyMesh>::FaceFace(m);
        vcg::tri::UpdateTopology<MyMesh>::VertexFace(m);
        MyMesh::VertexIterator vi;
        int checkflag;
        int i=0;
        int deg_vi;

        //index marking
        for(vi=m.vert.begin();vi!=m.vert.end();++vi){
            i++; (*vi).IMark() = i;
        }
        vi = m.vert.begin();
        int bi; 
        
        while(vi!=m.vert.end()){

            if((*vi).cFlags()) {vi++; continue;}
            //neighborhood check
            MyVertex* cvi;
            MyVertex* fNvi;
            MyVertex* Nvi;
            cvi = &(*vi);
            vcg::face::JumpingPos<MyFace> pos((*vi).VFp(),cvi);
            fNvi = pos.VFlip(); Nvi = fNvi;
            checkflag = 0;
            do {
                checkflag += (*Nvi).cFlags();
                pos.FlipF();
                pos.FlipE();
                Nvi = pos.VFlip();
            } while(Nvi != fNvi);
            if(checkflag) {vi++; continue;}
               
            // do extracting
            
            vcg::face::JumpingPos<MyFace> pos1((*vi).VFp(),cvi);
            fNvi = pos1.VFlip();
            Nvi = fNvi;
            CoordType pred(0,0,0);
            CoordType emb_vi(0,0,0);
            (*vi).Flags()=1;
            deg_vi = 0;
            do {
                deg_vi++;
                pred = pred + (*Nvi).cP();
                (*Nvi).Flags() = 1;
                pos1.FlipF();
                pos1.FlipE();
                Nvi = pos1.VFlip();
            } while(Nvi != fNvi);
            pred = pred/deg_vi;
            
            //std::cout << pred[1] << "  " << (*vi).P()[1] << std::endl;
            double di = CalculateDist(pred, (*vi).cP())/2;
            bi = (int)(di*std::pow(10,p_level)) % 2;

            std::cout << bi;
            double di_t = (di - (int)(di*std::pow(10,p_level)) / std::pow(10,p_level) - bi / (std::pow(10,p_level)) ) * 2;
            //double di_t = di + ((int)(di*std::pow(10,p_level)))/std::pow(10,p_level) + 1/(std::pow(10,p_level));
            emb_vi = pred + ((*vi).cP() - pred) * (di_t/(di));
            (*vi).P() = emb_vi;
            
            vi++;
            
        }
        std:: cout << std::endl;
        for(vi=m.vert.begin();vi!=m.vert.end();++vi){
            (*vi).Flags() = 0;       
        }
        std::cout << m.VN() << std::endl;
        //std::cout << i << std::endl;
        vcg::tri::io::Exporter<MyMesh>::Save(m, recov_mesh_name);
        if(ori_mesh_name){
            MyMesh ori;
            OpenMesh(ori, ori_mesh_name);
            PrintSNR(m, ori);
        }

    }
    else{
        std::cerr << "no options to be processed." << std::endl; 
    }

}

void OpenMesh( MyMesh &m, char *file_name )
{
    //std::cout << "1" <<std::endl;
    int err = vcg::tri::io::Importer<MyMesh>::Open(m, file_name);
    if(err){
        fprintf(stderr, "Error in loading %s: '%s'\n", file_name, vcg::tri::io::Importer<MyMesh>::ErrorMsg(err));
        if(vcg::tri::io::Importer<MyMesh>::ErrorCritical(err)) exit(-1);
    }
    // cleaning
}

int ParseArgs( int argc, char **argv ){
    //prg name
    argc--; argv++;
    //char *temp_addr = NULL;
    while(argc>0){
        if((*argv)[0] == '-'){ // option
            if(!strcmp(*argv, "-v")) {print_verbose = TRUE;}
            else if(!strcmp(*argv, "-ori")) {
                argv++; argc--;
                ori_mesh_name = *argv;
            }
            else if(!strcmp(*argv, "-stego")) {
                argv++; argc--;
                stego_mesh_name = *argv;
            }
            else if(!strcmp(*argv, "-recov")) {
                argv++; argc--;
                recov_mesh_name = *argv;
            }
            else if(!strcmp(*argv, "-embed")) {
                prc_embed = TRUE; prc_extract = FALSE;
            }
            else if(!strcmp(*argv, "-extract")) {
                prc_embed = FALSE; prc_extract = TRUE;
            }
            else if(!strcmp(*argv, "-wm")) {argv++;argc--; wm_file_name = *argv;}
            else if(!strcmp(*argv, "-m")) {argv++;argc--; p_level = atoi(*argv);}
            else{ fprintf(stderr, "Invalid program argument: %s\n", *argv); return 0; }
            argv++; argc--;
        }
        else { // mesh in-out name
            //if(!input_mesh_name) input_mesh_name = *argv;
            //else if(!output_mesh_name) output_mesh_name = *argv;
            //else 
            std::cerr << "Invalid program argument: " << *argv << std::endl;
            argv++; argc--; 
        }
    }
    if( prc_embed < 0 && prc_extract < 0 ){
        fprintf(stderr, "Usage: ./main <options> -embed or -extract \n\n");
        return FALSE;
    }
/*    if(!output_mesh_name){
        output_mesh_name = "temp.ply"; // check
    }

    QString file_name;
    QString short_file_name;
    QString out_directory;
    QString only_file_name;
    QDir dir(".");
    if( mode_embed ){
        QString file_name(QFileInfo(ori_mesh_name).absoluteFilePath());
        short_file_name = file_name.right(file_name.length() - file_name.lastIndexOf("/") - 1);
        only_file_name = short_file_name;
        short_file_name.truncate(short_file_name.lastIndexOf('.'));        
        out_directory = dir.absolutePath();
        //out_directory.truncate(file_name.lastIndexOf('/'));
        
        QString temp1 = QString("%1/%2/%3_m_%4_wm_%5.%6").arg(out_directory).arg("output/emb")
                                                        .arg(short_file_name).arg(m).arg(wm_file_name).arg("ply");
        QString temp2 = QString("%1/%2/%3_m_%4_wm_%5.%6").arg(out_directory).arg("output/log")
                                                        .arg(short_file_name).arg(m).arg(wm_file_name).arg("log");
        stego_mesh_name = temp1.toStdString().c_str();
        log_file_name = temp2.toStdString().c_str();
        
    }
    else { // mode_extract

    }
*/
    return TRUE;
}

int main( int argc, char** argv )
{
    //MyMesh                  M;
    unsigned long           elapsed_time; 
    
    int t0 = clock();
    
    if(!ParseArgs(argc,argv)) exit(-1);

    StartProcess();
    // OpenMesh(M);
    
    // //vcg::tri::io::ExporterPLY<MyMesh>::Save(M, "temp1.ply");
    // if(print_verbose){
    //     std::cout << "Mesh info: " << std::endl;
    //     fprintf(stdout, "  Input: '%s'\n\tvertices  %7i\n\tfaces    %7i\n", ori_mesh_name, M.vn, M.fn);
    // }
    // StartProcess(M);  
    //vcg::tri::io::ExporterPLY<MyMesh>::Save(M, "temp.ply");
    elapsed_time = clock() - t0;
    fprintf(stdout, "   Computation time  : %d ms\n", (int)(1000.0*elapsed_time/CLOCKS_PER_SEC));

    return 0;
}