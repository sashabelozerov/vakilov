/* 
 * File:   csvwriter.h
 * Author: wert
 *
 * Created on June 7, 2012, 12:38 AM
 */

#pragma once
#include <string>
#include <fstream>


class CSVWriter {
public:
	CSVWriter(std::string filename) {
		outFile.open(filename.c_str());
		if(!outFile.is_open())
			std::cerr << "unable to open file: " << filename << std::endl;
		
		outFile << "T"			<< "\t" 
				<< "eAvg"		<< "\t" 
				<< "eAvgError"	<< "\t"
				<< "e2Avg"		<< "\t"
				<< "e2AvgError"	<< "\t"
				<< "mAvg"		<< "\t" 
				<< "mAvgError"	<< "\t"
				<< "m2Avg"		<< "\t"
				<< "m2AvgError"	<< "\t"
				<< "c"			<< "\t"
				<< "cError"		<< "\t"
				<< "x"			<< "\t"
				<< "xError"		<< std::endl;
	}
	
	virtual ~CSVWriter(){
		if(outFile.is_open())
			outFile.close();
	}
	
	void write(	double T, 
				double eAvg, 
				double eAvgError, 
				double e2Avg,
				double e2AvgError,
				double mAvg, 
				double mAvgError,
				double m2Avg, 
				double m2AvgError,
				double c,
				double cError,
				double x,
				double xError) {

		outFile << T			<< "\t" 
				<< eAvg			<< "\t" 
				<< eAvgError	<< "\t"
				<< e2Avg		<< "\t"
				<< e2AvgError	<< "\t"
				<< mAvg			<< "\t" 
				<< mAvgError	<< "\t"
				<< m2Avg		<< "\t"
				<< m2AvgError	<< "\t"
				<< c			<< "\t"
				<< cError		<< "\t"
				<< x			<< "\t"
				<< xError		<< std::endl;
	}
	
private:
	std::string filename;
	std::ofstream outFile;
};