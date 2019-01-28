#ifndef UtilsH
#define UtilsH

#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include <cstring>

class Utils{
  public:
    static std::vector<double> getValuesFromCsv(std::string filename, std::string start_date, std::string end_date,int value_column){
      std::ifstream file = std::ifstream(filename,std::ifstream::in);
  		if( file.fail() ) std::cout << "Couldn't open the file";
      std::vector<double> values;
      bool in_range=false;
  		while(!file.eof()){
          std::string line;
  		    std::getline(file,line);
          if(line=="") continue;
          std::vector<std::string> tokens=parseLine(line.c_str(),",");
          std::string date=tokens[0];
          if(date==start_date) in_range=true;
          if(date==end_date) in_range=false;
          if(in_range){
            values.push_back(std::stod(tokens[value_column]));
          }
  		}
      file.close();

      return values;
    }


    static std::vector<std::string> parseLine(const char* constLine,const char* delimiter){
        char* line=strdup(constLine);
        std::vector<std::string> tokens;
        char* token=strtok(line,delimiter);
        tokens.push_back(token);
        while( (token=strtok(NULL,delimiter)) ){
            tokens.push_back(token);
        }

        return tokens;
    }

    static std::vector<double> getAverageTemperaturesFromCsv(std::string filename, std::string start_date, std::string end_date){///#in Kelvin#TODO:change name to meanTemperature
        std::vector<double> temperatures=getValuesFromCsv(filename,start_date,end_date,2);
        for(unsigned int i=0;i<temperatures.size();i++) temperatures[i]+=273.15;
        return temperatures;
    }

    static std::vector<double> getPrecipitationsFromCsv(std::string filename, std::string start_date, std::string end_date){//in mm
        return getValuesFromCsv(filename,start_date,end_date,4);
    }

    static std::vector<double> getRelativeHumidityFromCsv(std::string filename, std::string start_date, std::string end_date){//in percentage
        return getValuesFromCsv(filename,start_date,end_date,5);
    }

    static unsigned int getDaysFromCsv(std::string filename, std::string start_date, std::string end_date){//convenience method to get amount of days between two dates
        return getValuesFromCsv(filename,start_date,end_date,0).size();
    }
};

#endif
