#include <crsRead/MCorsikaReader.h>
#include <crs/TSubBlock.h>
#include <crs/MRunHeader.h>
#include <crs/MEventHeader.h>
#include <crs/MEventEnd.h>
#include <crs/MParticleBlock.h>
#include <crs/MLongitudinalBlock.h>
#include <crs/MParticle.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <fstream>  // Adicionado para salvar dados em arquivos

#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <unistd.h>

std::string usage();
int getOptions(int argc, char** argv, int& n, std::string& outFilename);

int main( int argc, char **argv )
{
  int shIndex = 1;
  std::string outFilename = "footprint_data.csv";  // Nome padrão do arquivo de saída
  
  const int nOptions = getOptions(argc, argv, shIndex, outFilename);  // Passa o nome do arquivo
  if(nOptions < 0 || argc-nOptions < 1 )
  {
    std::cout << usage() << std::endl;
    return 1;
  }

  /*
   * -----------------------------------------------
   * Reading the file
   * -----------------------------------------------
   */

  std::string fname(argv[nOptions]);

  crsRead::MCorsikaReader cr(fname, 1);
  crs::MRunHeader run;

  Double_t hs = 2000.0;
  Int_t nBins = hs;

  TH2D* xy_tot = new TH2D("xy_tot","xy_tot",nBins,-hs,hs,nBins,-hs,hs);
  TH2D* xy_e = new TH2D("xy_e","xy_e",nBins,-hs,hs,nBins,-hs,hs);
  TH2D* xy_c = new TH2D("xy_c","xy_c",nBins,-hs,hs,nBins,-hs,hs);
  TH2D* xy_m = new TH2D("xy_m","xy_m",nBins,-hs,hs,nBins,-hs,hs);
  TH2D* xy_n = new TH2D("xy_n","xy_n",nBins,-hs,hs,nBins,-hs,hs);
  
  Int_t i=0;
  Double_t logE = 0.0;
  Double_t obsLev = 0.0;

  std::ofstream outFile(outFilename);  // Usar o nome do arquivo fornecido
  outFile << "X,Y,Type\n";  // Cabeçalho para o arquivo CSV
  
  while( cr.GetRun( run ) ) // Loop over the runs. Usually one
  {
    crs::MEventHeader shower;
    while( cr.GetShower( shower ) )
    {
      i++;
      if(i != shIndex)
        continue;

      std::cout << "-----------------------------" << std::endl;
      std::cout << "plot of the shower nr: " << i << std::endl;
      logE = log10( shower.GetEnergy() );
      obsLev = shower.GetObservationHeight( shower.GetNObservationLevels()-1 )/100.0;
      std::cout << "shower obs level: " << obsLev << " m" << std::endl;
      
      crs::TSubBlock sub_block;
      while( cr.GetData( sub_block ) )
      {
        crs::MParticleBlock data(sub_block);
        if( data.GetBlockType() == crs::TSubBlock::ePARTDATA )
        {
          for( crs::MParticleBlock::ParticleListConstIterator p_it = data.FirstParticle();
               p_it != data.LastParticle();
               ++p_it )
          {
            crs::MParticle part(*p_it);

            Int_t pID = part.GetParticleID();

            // Salvando as partículas em um arquivo CSV com base no tipo
            if( pID < 4 )
            {
              xy_e->Fill(part.GetX()/100.0, part.GetY()/100.0);
              outFile << part.GetX()/100.0 << "," << part.GetY()/100.0 << ",electron\n";  // Salva elétron
            }

            if( pID == 5 || pID == 6 )
            {
              xy_m->Fill(part.GetX()/100.0, part.GetY()/100.0);
              outFile << part.GetX()/100.0 << "," << part.GetY()/100.0 << ",muon\n";  // Salva múon
            }

            if( part.IsCherenkov() )
            {
              xy_c->Fill(part.GetX()/100.0, part.GetY()/100.0);
              outFile << part.GetX()/100.0 << "," << part.GetY()/100.0 << ",cherenkov\n";  // Salva fótons de Cherenkov
            }

            if( part.IsNucleus() || ( pID >6 && pID <200 ) && !part.IsMuonProductionInfo())
            {
              xy_n->Fill(part.GetX()/100.0, part.GetY()/100.0);
              outFile << part.GetX()/100.0 << "," << part.GetY()/100.0 << ",nucleus\n";  // Salva núcleos
            }

            if( !part.IsMuonProductionInfo() )
            {
              xy_tot->Fill(part.GetX()/100.0, part.GetY()/100.0);
              outFile << part.GetX()/100.0 << "," << part.GetY()/100.0 << ",total\n";  // Salva todas as partículas
            }
          }
        }
      }
      break;
    }
  }

  outFile.close();  // Fecha o arquivo CSV ao final da execução

  if ( i < shIndex )
  {
    std::cout << "ERROR: The file contains " << i << " showers" << std::endl;
    return -1;
  }

  std::cout << "-----------------------------" << std::endl;

  return 0;
}

std::string usage()
{
  std::string usg = "usage: ./energySpectra <corsika file name> -o <output csv filename>";
  return usg;
}

int getOptions( int argc, char** argv, int &n, std::string& outFilename )
{
  int c;
  while ( (c = getopt(argc, argv, "n:o:h")) != -1 )
  {
    switch(c)
    {
    case 'n':
      n = std::atoi(optarg);
      break;

    case 'o':  // Novo argumento para o nome do arquivo de saída
      outFilename = optarg;
      break;
      
    case 'h':
      return -2;

    default:
      return -2;
    }
  }
  return optind;
}

