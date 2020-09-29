//--------------------------------------------------------------------------
// This file is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This file is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with This file. If not, see <http://www.gnu.org/licenses/>.
//--------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <complex>
#include <iomanip>
#include <assert.h>
#include <algorithm>
#include <stdio.h>
#include <string.h>

#define _USE_MATH_DEFINES
#include <math.h>

enum StationType{
	MT = 0,
	VTF,
	PT,
	HTF,
	NMT,
	NMT2,
	APP_RES_AND_PHS,
	NMT2_APP_RES_AND_PHS,
	ADDITINAL_OUTPUT_POINT
};

enum ImpedanceTensorComponentIndex{
	XX = 0,
	XY,
	YX,
	YY,
};

struct ImpedanceTensor{
	std::complex<double> Z[4];
};

struct VerticalMagneticTransferFunction{
	std::complex<double> TZ[2];
};

struct ApparentResistivityAndPhase{
	double apparentResistivity[4];
	double phase[4];
};

struct DistortionMatrix{
	double C[4];
};

struct MTData{
	double freq;
	ImpedanceTensor Cal;
	ImpedanceTensor Res;
	ImpedanceTensor Obs;
	ImpedanceTensor Err;
};

struct VTFData{
	double freq;
	VerticalMagneticTransferFunction Cal;
	VerticalMagneticTransferFunction Res;
	VerticalMagneticTransferFunction Obs;
	VerticalMagneticTransferFunction Err;
};

struct ApparentResistivityAndPhaseData{
	double freq;
	ApparentResistivityAndPhase Cal;
	ApparentResistivityAndPhase Res;
	ApparentResistivityAndPhase Obs;
	ApparentResistivityAndPhase Err;
};

struct MTTrueError{
	double freq;
	std::pair<double, double> error[4];
};

struct VTFTrueError{
	double freq;
	std::pair<double, double> error[2];
};

struct NumDataAndSumResidual{
	int numData;
	double sumResidual;
};

// Type of galvanic distortion
enum TypeOfDistortion{
	NO_DISTORTION = 0,
	ESTIMATE_DISTORTION_MATRIX_DIFFERENCE,
	ESTIMATE_GAINS_AND_ROTATIONS,
	ESTIMATE_GAINS_ONLY,
};

const double MU0 = 4.0 * M_PI * 1.0e-7;
const double RAD2DEG = 180.0 / M_PI;
const double DEG2RAD =  M_PI /180.0;
const double FACTOR = 1.0;

const std::string m_componentIndex[6] = { "xx", "xy", "yx", "yy", "zx", "zy" };
int m_stationTypeCur(-1);
bool m_readTrueErrorFile(false);
bool m_outputCSV(false);
bool m_isImpedanceTensorConvertedToAppResAndPhase(false);
int m_typeOfDistortion(NO_DISTORTION);
std::vector< std::pair<int, MTData> > m_MTDataListAll;
std::vector< std::pair<int, VTFData> > m_VTFDataListAll;
std::vector< std::pair<int, MTData> > m_NMT2DataListAll;
std::vector< std::pair<int, ApparentResistivityAndPhaseData> > m_AppResAndPhsDataListAll;
std::vector< std::pair<int, ApparentResistivityAndPhaseData> > m_NMT2AppResAndPhsDataListAll;
std::vector< std::pair<int, MTTrueError> > m_MTTrueErrorListAll;
std::vector< std::pair<int, VTFTrueError> > m_VTFTrueErrorListAll;
std::vector< std::pair<int, MTTrueError> > m_NMT2TrueErrorListAll;
std::vector< std::pair<int, MTTrueError> > m_AppResAndPhsTrueErrorListAll;
std::vector< std::pair<int, MTTrueError> > m_NMT2AppResAndPhsTrueErrorListAll;
std::map< int, StationType > m_SiteIDToSiteType;
std::map< int, DistortionMatrix > m_distortionMatrixList;
std::map< int, std::string > m_SiteIDToSiteName;

void readResult( const int iterationNumber, const int numPE );
void readControlData( const std::string& fileName );
void readTrueError( const std::string& trueErrorFileName );
void readDistortionMatrix( const int iterationNumber );
void readRelationSiteIDBetweenSiteName( const std::string& relationFile );
void calcTrueRMS( const std::string& trueErrorFileName );
std::string convertSiteIDToSiteName( const int siteID );
void writeResult();
void writeResultVTF();
void writeResultMT();
void writeResultNMT2();
void writeResultAppResAndPhs();
void writeResultNMT2AppResAndPhs();
bool pairCompareMTData( const std::pair<int, MTData>& left, const std::pair<int, MTData>& right );
bool pairCompareVTFData( const std::pair<int, VTFData>& left, const std::pair<int, VTFData>& right );
bool pairCompareAppResAndPhsData( const std::pair<int, ApparentResistivityAndPhaseData>& left, const std::pair<int, ApparentResistivityAndPhaseData>& right );
void calcUndistortedImpedanceTensor( const int siteID, const ImpedanceTensor& Z, ImpedanceTensor& ZWithoutDistortion);
void calcUndistortedApparentResistiivtyAndPhase( const int siteID, const double freq, const ApparentResistivityAndPhase& appResAndPhs, ApparentResistivityAndPhase& appResAndPhsWithoutDistortion );
int getSiteTypeFromSiteID( const int siteID );
void addNumDataAndSumResidual( const int siteID, std::map<int, NumDataAndSumResidual> & numDataAndSumResidual, const double residual );
std::vector< std::pair<int, MTTrueError> >::const_iterator getIteratorToMTTrueError( const int siteID, const double freq );
std::vector< std::pair<int, VTFTrueError> >::const_iterator getIteratorToVTFTrueError( const int siteID, const double freq );
std::vector< std::pair<int, MTTrueError> >::const_iterator getIteratorToNMT2TrueError( const int siteID, const double freq );
std::vector< std::pair<int, MTTrueError> >::const_iterator getIteratorToAppResAndPhsTrueError( const int siteID, const double freq );
std::vector< std::pair<int, MTTrueError> >::const_iterator getIteratorToNMT2AppResAndPhsTrueError( const int siteID, const double freq );
void isSameError( const int siteID, const double freq, const int index, const double errReal, const double errImag );
void writeToOfstream( std::ofstream& ofs, const int width, const std::string& str, const bool isCSV );
void writeToOfstream( std::ofstream& ofs, const int width, const int value, const bool isCSV );
void writeToOfstream( std::ofstream& ofs, const int width, const double value, const bool isCSV );

int main( int argc, char* argv[] ){
	if( argc < 3 ){
		std::cerr << "You must specify iteration number and process number !!" << std::endl;
		exit(1);
	}
	readResult( atoi(argv[1]), atoi(argv[2]) );
	for( int i = 3; i < argc; ++i ){
		if( strcmp(argv[i], "-name") == 0 ) {
			readRelationSiteIDBetweenSiteName( argv[i+1] );
		}
		else if( strcmp(argv[i], "-err") == 0 ) {
			calcTrueRMS( argv[i+1] );
			m_readTrueErrorFile = true;
		}
		else if( strcmp(argv[i], "-csv") == 0 ) {
			m_outputCSV = true;
		}
		else if( strcmp(argv[i], "-undist") == 0 ) {
			readControlData( "control.dat" );
			readDistortionMatrix( atoi(argv[1]) );
		}
		else if( strcmp(argv[i], "-appphs") == 0 ) {
			std::cout << "Impedance tensors are converted to apparent resistivity and phase." << std::endl;
			m_isImpedanceTensorConvertedToAppResAndPhase = true;
		}
	}
	writeResult();
	return 0;
}

void readResult( const int iterationNumber, const int numPE ){

	for( int pe = 0; pe < numPE; ++pe )
	{
		std::ostringstream resultFileName;
		resultFileName << "result_" << pe << "_iter" << iterationNumber << ".csv";
		std::ifstream ifs( resultFileName.str().c_str(), std::ios::in );
		if( ifs.fail() )
		{
			std::cerr << "File open error : " << resultFileName.str() << " !!" << std::endl;
			exit(1);
		}
		std::cout << "Read result from " << resultFileName.str() << std::endl;
		
		std::string str;
		while(getline(ifs, str))
		{
			std::istringstream ss(str);
			std::string sbuf;
			getline(ss, sbuf, ',');// Read first data of this line
			if( sbuf == "MT" ){
				m_stationTypeCur = MT;
				continue;
			}else if(  sbuf == "VTF"  ){
				m_stationTypeCur = VTF;
				continue;
			}else if(  sbuf == "PT"  ){
				m_stationTypeCur = PT;
				continue;
			}else if(  sbuf == "HTF"  ){
				m_stationTypeCur = HTF;
				continue;
			}else if(  sbuf == "NMT"  ){
				m_stationTypeCur = NMT;
				continue;
			}else if(  sbuf == "NMT2"  ){
				m_stationTypeCur = NMT2;
				continue;
			}else if(  sbuf == "APP_RES_AND_PHS"  ){
				m_stationTypeCur = APP_RES_AND_PHS;
				continue;
			}else if(  sbuf == "NMT2_APP_RES_AND_PHS"  ){
				m_stationTypeCur = NMT2_APP_RES_AND_PHS;
				continue;
			}else if(  sbuf == "ADDITINAL_OUTPUT_POINT"  ){
				m_stationTypeCur = ADDITINAL_OUTPUT_POINT;
				continue;
			}else if(  sbuf.find("StaID") != std::string::npos ){
				continue;
			}
			// Read response functions
			//const int statID = std::stoi(sbuf);
			std::istringstream ssStatID(sbuf);
			int statID(-1);
			ssStatID >> statID;
			char cbuf;
			if( m_stationTypeCur == MT ) {
				MTData data;
				ss >> data.freq >> cbuf;
				double dbufReal(0.0);
				double dbufImag(0.0);
				for( int i = 0; i < 4; ++i ){
					ss >> dbufReal >> cbuf >> dbufImag >> cbuf;
					data.Cal.Z[i] = std::complex<double>(dbufReal, dbufImag);
				}
				for( int i = 0; i < 4; ++i ){
					ss >> dbufReal >> cbuf >> dbufImag >> cbuf;
					data.Res.Z[i] = std::complex<double>(dbufReal, dbufImag);
				}
				for( int i = 0; i < 4; ++i ){
					ss >> dbufReal >> cbuf >> dbufImag >> cbuf;
					data.Obs.Z[i] = std::complex<double>(dbufReal, dbufImag);
				}
				for( int i = 0; i < 4; ++i ){
					ss >> dbufReal >> cbuf >> dbufImag >> cbuf;
					data.Err.Z[i] = std::complex<double>(dbufReal, dbufImag);
				}
				m_MTDataListAll.push_back( std::make_pair(statID, data) );
				m_SiteIDToSiteType.insert( std::make_pair(statID, MT) );
			}else if( m_stationTypeCur == VTF ) {
				VTFData data;
				ss >> data.freq >> cbuf;
				double dbufReal(0.0);
				double dbufImag(0.0);
				for( int i = 0; i < 2; ++i ){
					ss >> dbufReal >> cbuf >> dbufImag >> cbuf;
					data.Cal.TZ[i] = std::complex<double>(dbufReal, dbufImag);
				}
				for( int i = 0; i < 2; ++i ){
					ss >> dbufReal >> cbuf >> dbufImag >> cbuf;
					data.Res.TZ[i] = std::complex<double>(dbufReal, dbufImag);
				}
				for( int i = 0; i < 2; ++i ){
					ss >> dbufReal >> cbuf >> dbufImag >> cbuf;
					data.Obs.TZ[i] = std::complex<double>(dbufReal, dbufImag);
				}
				for( int i = 0; i < 2; ++i ){
					ss >> dbufReal >> cbuf >> dbufImag >> cbuf;
					data.Err.TZ[i] = std::complex<double>(dbufReal, dbufImag);
				}
				m_VTFDataListAll.push_back( std::make_pair(statID, data) );
				m_SiteIDToSiteType.insert( std::make_pair(statID, VTF ) );
			}else if( m_stationTypeCur == NMT2 ) {
				MTData data;
				ss >> data.freq >> cbuf;
				double dbufReal(0.0);
				double dbufImag(0.0);
				for( int i = 0; i < 4; ++i ){
					ss >> dbufReal >> cbuf >> dbufImag >> cbuf;
					data.Cal.Z[i] = std::complex<double>(dbufReal, dbufImag);
				}
				for( int i = 0; i < 4; ++i ){
					ss >> dbufReal >> cbuf >> dbufImag >> cbuf;
					data.Res.Z[i] = std::complex<double>(dbufReal, dbufImag);
				}
				for( int i = 0; i < 4; ++i ){
					ss >> dbufReal >> cbuf >> dbufImag >> cbuf;
					data.Obs.Z[i] = std::complex<double>(dbufReal, dbufImag);
				}
				for( int i = 0; i < 4; ++i ){
					ss >> dbufReal >> cbuf >> dbufImag >> cbuf;
					data.Err.Z[i] = std::complex<double>(dbufReal, dbufImag);
				}
				m_NMT2DataListAll.push_back( std::make_pair(statID, data) );
				m_SiteIDToSiteType.insert( std::make_pair(statID, NMT2) );
			}else if( m_stationTypeCur == APP_RES_AND_PHS ) {
				ApparentResistivityAndPhaseData data;
				ss >> data.freq >> cbuf;
				double dbufAppRes(0.0);
				double dbufPhs(0.0);
				for( int i = 0; i < 4; ++i ){
					ss >> dbufAppRes >> cbuf >> dbufPhs >> cbuf;
					data.Cal.apparentResistivity[i] = dbufAppRes;
					data.Cal.phase[i] = dbufPhs;
				}
				for( int i = 0; i < 4; ++i ){
					ss >> dbufAppRes >> cbuf >> dbufPhs >> cbuf;
					data.Res.apparentResistivity[i] = dbufAppRes;
					data.Res.phase[i] = dbufPhs;
				}
				for( int i = 0; i < 4; ++i ){
					ss >> dbufAppRes >> cbuf >> dbufPhs >> cbuf;
					data.Obs.apparentResistivity[i] = dbufAppRes;
					data.Obs.phase[i] = dbufPhs;
				}
				for( int i = 0; i < 4; ++i ){
					ss >> dbufAppRes >> cbuf >> dbufPhs >> cbuf;
					data.Err.apparentResistivity[i] = dbufAppRes;
					data.Err.phase[i] = dbufPhs;
				}
				m_AppResAndPhsDataListAll.push_back( std::make_pair(statID, data) );
				m_SiteIDToSiteType.insert( std::make_pair(statID, APP_RES_AND_PHS) );
			}else if( m_stationTypeCur == NMT2_APP_RES_AND_PHS ) {
				ApparentResistivityAndPhaseData data;
				ss >> data.freq >> cbuf;
				double dbufAppRes(0.0);
				double dbufPhs(0.0);
				for( int i = 0; i < 4; ++i ){
					ss >> dbufAppRes >> cbuf >> dbufPhs >> cbuf;
					data.Cal.apparentResistivity[i] = dbufAppRes;
					data.Cal.phase[i] = dbufPhs;
				}
				for( int i = 0; i < 4; ++i ){
					ss >> dbufAppRes >> cbuf >> dbufPhs >> cbuf;
					data.Res.apparentResistivity[i] = dbufAppRes;
					data.Res.phase[i] = dbufPhs;
				}
				for( int i = 0; i < 4; ++i ){
					ss >> dbufAppRes >> cbuf >> dbufPhs >> cbuf;
					data.Obs.apparentResistivity[i] = dbufAppRes;
					data.Obs.phase[i] = dbufPhs;
				}
				for( int i = 0; i < 4; ++i ){
					ss >> dbufAppRes >> cbuf >> dbufPhs >> cbuf;
					data.Err.apparentResistivity[i] = dbufAppRes;
					data.Err.phase[i] = dbufPhs;
				}
				m_NMT2AppResAndPhsDataListAll.push_back( std::make_pair(statID, data) );
				m_SiteIDToSiteType.insert( std::make_pair(statID, NMT2_APP_RES_AND_PHS) );
			}else if( m_stationTypeCur == ADDITINAL_OUTPUT_POINT ) {
				break;
			}else{
				std::cerr << "Unsupported station type : " << m_stationTypeCur << std::endl; 
				break;
			}
		}
	}
}

void readControlData( const std::string& fileName ){

	std::ifstream ifs( fileName.c_str() );
	if( ifs.fail() )
	{
		std::cerr << "File open error : " << fileName << " !!" << std::endl;
		exit(1);
	}

	bool found(false);
	while(!ifs.eof())
	{
		std::string line;
		ifs >> line;

		double dbuf(0.0);
		int ibuf(0);
		if( line.substr(0,10).compare("DISTORTION") == 0 ){
			ifs >> ibuf;

			if( ibuf != NO_DISTORTION && 
				ibuf != ESTIMATE_DISTORTION_MATRIX_DIFFERENCE &&
				ibuf != ESTIMATE_GAINS_AND_ROTATIONS &&
				ibuf != ESTIMATE_GAINS_ONLY ){
				std::cerr << "Error : Wrong type ID is written below DISTORTION : " << ibuf << std::endl;
				exit(1);
				break;
			}

			m_typeOfDistortion = ibuf;

			found = true;
		}else if( line.substr(0,3).compare("END") == 0 ){
			break;
		}		
	}
	
	if( found ){
		std::cout << "Type of distortion : " << m_typeOfDistortion << std::endl;
	}else{
		std::cerr << "DISTORTION is not found in " << fileName << " !!" << std::endl;
	}

	ifs.close();
}

void readTrueError( const std::string& trueErrorFileName ){

	std::ifstream ifs( trueErrorFileName.c_str() );
	if( ifs.fail() )
	{
		std::cerr << "File open error : " << trueErrorFileName << " !!" << std::endl;
		exit(1);
	}

	std::string str;
	while(getline(ifs, str))
	{
		std::istringstream ss(str);
#ifdef _DEBUG
		std::cout << str << std::endl;
#endif
		int statID(-1);
		ss >> statID;
		if( statID < 0 ){
			break;
		}
		const int type = getSiteTypeFromSiteID(statID);
		if( type == MT ) {
			MTTrueError trueError;
			ss >> trueError.freq;
			double dbufReal(0.0);
			double dbufImag(0.0);
			for( int i = 0; i < 4; ++i ){
				ss >> dbufReal >> dbufImag;
				trueError.error[i] = std::make_pair(dbufReal, dbufImag);
			}
			m_MTTrueErrorListAll.push_back( std::make_pair(statID, trueError) );
		}else if( type == VTF ) {
			VTFTrueError trueError;
			ss >> trueError.freq;
			double dbufReal(0.0);
			double dbufImag(0.0);
			for( int i = 0; i < 2; ++i ){
				ss >> dbufReal >> dbufImag;
				trueError.error[i] = std::make_pair(dbufReal, dbufImag);
			}
			m_VTFTrueErrorListAll.push_back( std::make_pair(statID, trueError) );
		}else if( type == NMT2 ) {
			MTTrueError trueError;
			ss >> trueError.freq;
			double dbufReal(0.0);
			double dbufImag(0.0);
			for( int i = 0; i < 4; ++i ){
				ss >> dbufReal >> dbufImag;
				trueError.error[i] = std::make_pair(dbufReal, dbufImag);
			}
			m_NMT2TrueErrorListAll.push_back( std::make_pair(statID, trueError) );
		}else if( type == APP_RES_AND_PHS ) {
			MTTrueError trueError;
			ss >> trueError.freq;
			double dbufAppRes(0.0);
			double dbufPhs(0.0);
			for( int i = 0; i < 4; ++i ){
				ss >> dbufAppRes >> dbufPhs;
				trueError.error[i] = std::make_pair(dbufAppRes, dbufPhs);
			}
			m_AppResAndPhsTrueErrorListAll.push_back( std::make_pair(statID, trueError) );
		}else if( type == NMT2_APP_RES_AND_PHS ) {
			MTTrueError trueError;
			ss >> trueError.freq;
			double dbufAppRes(0.0);
			double dbufPhs(0.0);
			for( int i = 0; i < 4; ++i ){
				ss >> dbufAppRes >> dbufPhs;
				trueError.error[i] = std::make_pair(dbufAppRes, dbufPhs);
			}
			m_NMT2AppResAndPhsTrueErrorListAll.push_back( std::make_pair(statID, trueError) );
		}else{
			std::cerr << "Unsupported station type : " << type << std::endl; 
			break;
		}
	}
	
	ifs.close();

}

void readDistortionMatrix( const int iterationNumber ){

	if( m_typeOfDistortion == NO_DISTORTION ){
		return;
	}

	std::ostringstream oss;
	oss << "distortion_iter" << iterationNumber << ".dat";
	const std::string fileName = oss.str();
	std::ifstream ifs( fileName.c_str() );
	if( ifs.fail() )
	{
		std::cout << fileName << " is not found." << std::endl;
		exit(1);
	}
	std::cout << "Read distortion matrix from " << fileName << std::endl;
	
	int numSite(-1);
	ifs >> numSite;

	for( int isite = 0; isite < numSite; ++isite ){
		int siteID(-1);
		ifs >> siteID;
		DistortionMatrix distorionMatrix;
		if( m_typeOfDistortion == ESTIMATE_DISTORTION_MATRIX_DIFFERENCE ){
			for( int i = 0; i < 4; ++i ){
				ifs >> distorionMatrix.C[i];
			}
			distorionMatrix.C[XX] += 1.0;
			distorionMatrix.C[YY] += 1.0;
		}
		else if( m_typeOfDistortion == ESTIMATE_GAINS_AND_ROTATIONS ){
			double dbuf[4];
			for( int i = 0; i < 4; ++i ){
				ifs >> dbuf[i];
			}		
			const double gX = pow( 10.0, dbuf[0] );
			const double gY = pow( 10.0, dbuf[1] );
			const double betaX = dbuf[2] * DEG2RAD;
			const double betaY = dbuf[3] * DEG2RAD;
			distorionMatrix.C[XX] =   gX * cos( betaX );
			distorionMatrix.C[XY] = - gY * sin( betaY );
			distorionMatrix.C[YX] =   gX * sin( betaX );
			distorionMatrix.C[YY] =   gY * cos( betaY );
		}
		else if( m_typeOfDistortion == ESTIMATE_GAINS_ONLY ){
			double dbuf[2];
			for( int i = 0; i < 2; ++i ){
				ifs >> dbuf[i];
			}		
			const double gX = pow( 10.0, dbuf[0] );
			const double gY = pow( 10.0, dbuf[1] );
			distorionMatrix.C[XX] =  gX;
			distorionMatrix.C[XY] = 0.0;
			distorionMatrix.C[YX] = 0.0;
			distorionMatrix.C[YY] =  gY;
		}
		else{
			std::cerr << "Unsupported type of distortion : " << m_typeOfDistortion << std::endl;
		}
		m_distortionMatrixList.insert( std::make_pair( siteID, distorionMatrix ) );
		int ibuf(-1);
		ifs >> ibuf;
	}

	ifs.close();

}

void readRelationSiteIDBetweenSiteName( const std::string& relationFile ){
	
	std::ifstream ifs( relationFile.c_str() );
	if( ifs.fail() )
	{
		std::cerr << "File open error : " << relationFile << " !!" << std::endl;
		exit(1);
	}
	std::cout << "Read relation between site ID and site name from " << relationFile << std::endl;

	std::string str;
	while(getline(ifs, str))
	{
		std::istringstream ss(str);
		int id(-1);
		std::string name;
		ss >> id >> name;
		m_SiteIDToSiteName.insert( std::make_pair( id, name ) );
	}

	ifs.close();

}

void calcTrueRMS( const std::string& trueErrorFileName ){

	readTrueError( trueErrorFileName );

	std::map<int, NumDataAndSumResidual> numDataAndSumResidualEachSite;

	// MT
	for( std::vector<  std::pair<int, MTData> >::const_iterator itrData = m_MTDataListAll.begin(); itrData != m_MTDataListAll.end(); ++itrData ){
		const int siteID = itrData->first;
		const double freq = itrData->second.freq;
		std::vector< std::pair<int, MTTrueError> >::const_iterator itrErr = getIteratorToMTTrueError( siteID, freq );
		for( int i = 0; i < 4; ++i ){
			// Real part
			if( itrErr->second.error[i].first > 0.0 ){
				const double residual = ( itrData->second.Cal.Z[i].real() - itrData->second.Obs.Z[i].real() ) / itrErr->second.error[i].first;
				addNumDataAndSumResidual(siteID, numDataAndSumResidualEachSite, residual);
			}
			// Imaginary part
			if( itrErr->second.error[i].second > 0.0 ){
				const double residual = ( itrData->second.Cal.Z[i].imag() - itrData->second.Obs.Z[i].imag() ) / itrErr->second.error[i].second;
				addNumDataAndSumResidual(siteID, numDataAndSumResidualEachSite, residual);
			}
		}
	}
	
	// VTF
	// Not checked
	for( std::vector<  std::pair<int, VTFData> >::const_iterator itrData = m_VTFDataListAll.begin(); itrData != m_VTFDataListAll.end(); ++itrData ){
		const int siteID = itrData->first;
		const double freq = itrData->second.freq;
		std::vector< std::pair<int, VTFTrueError> >::const_iterator itrErr = getIteratorToVTFTrueError( siteID, freq );
		for( int i = 0; i < 2; ++i ){
			// Real part
			if( itrErr->second.error[i].first > 0.0 ){
				const double residual = ( itrData->second.Cal.TZ[i].real() - itrData->second.Obs.TZ[i].real() ) / itrErr->second.error[i].first;
				addNumDataAndSumResidual(siteID, numDataAndSumResidualEachSite, residual);
			}
			// Imaginary part
			if( itrErr->second.error[i].second > 0.0 ){
				const double residual = ( itrData->second.Cal.TZ[i].imag() - itrData->second.Obs.TZ[i].imag() ) / itrErr->second.error[i].second;
				addNumDataAndSumResidual(siteID, numDataAndSumResidualEachSite, residual);
			}
		}
	}

	// NMT2
	for( std::vector<  std::pair<int, MTData> >::const_iterator itrData = m_NMT2DataListAll.begin(); itrData != m_NMT2DataListAll.end(); ++itrData ){
		const int siteID = itrData->first;
		const double freq = itrData->second.freq;
		std::vector< std::pair<int, MTTrueError> >::const_iterator itrErr = getIteratorToNMT2TrueError( siteID, freq );
		for( int i = 0; i < 4; ++i ){
			// Real part
			if( itrErr->second.error[i].first > 0.0 ){
				const double residual = ( itrData->second.Cal.Z[i].real() - itrData->second.Obs.Z[i].real() ) / itrErr->second.error[i].first;
				addNumDataAndSumResidual(siteID, numDataAndSumResidualEachSite, residual);
			}
			// Imaginary part
			if( itrErr->second.error[i].second > 0.0 ){
				const double residual = ( itrData->second.Cal.Z[i].imag() - itrData->second.Obs.Z[i].imag() ) / itrErr->second.error[i].second;
				addNumDataAndSumResidual(siteID, numDataAndSumResidualEachSite, residual);
			}
		}
	}

	// APP_RES_AND_PHS
	for( std::vector<  std::pair<int, ApparentResistivityAndPhaseData> >::const_iterator itrData = m_AppResAndPhsDataListAll.begin();
		itrData != m_AppResAndPhsDataListAll.end(); ++itrData ){
		const int siteID = itrData->first;
		const double freq = itrData->second.freq;
		std::vector< std::pair<int, MTTrueError> >::const_iterator itrErr = getIteratorToAppResAndPhsTrueError( siteID, freq );
		for( int i = 0; i < 4; ++i ){
			// Apparent resistivity
			if( itrErr->second.error[i].first > 0.0 ){
				const double residual = log10( itrData->second.Cal.apparentResistivity[i] / itrData->second.Obs.apparentResistivity[i] )
					/ log10 ( 1.0 + itrErr->second.error[i].first / itrData->second.Obs.apparentResistivity[i] );
				addNumDataAndSumResidual(siteID, numDataAndSumResidualEachSite, residual);
			}
			// Phase
			if( itrErr->second.error[i].second > 0.0 ){
				const double residual = ( itrData->second.Cal.phase[i] - itrData->second.Obs.phase[i] ) / itrErr->second.error[i].second;
				addNumDataAndSumResidual(siteID, numDataAndSumResidualEachSite, residual);
			}
		}
	}

	// NMT2_APP_RES_AND_PHS
	for( std::vector<  std::pair<int, ApparentResistivityAndPhaseData> >::const_iterator itrData = m_NMT2AppResAndPhsDataListAll.begin();
		itrData != m_NMT2AppResAndPhsDataListAll.end(); ++itrData ){
		const int siteID = itrData->first;
		const double freq = itrData->second.freq;
		std::vector< std::pair<int, MTTrueError> >::const_iterator itrErr = getIteratorToNMT2AppResAndPhsTrueError( siteID, freq );
		for( int i = 0; i < 4; ++i ){
			// Apparent resistivity
			if( itrErr->second.error[i].first > 0.0 ){
				const double residual = log10( itrData->second.Cal.apparentResistivity[i] / itrData->second.Obs.apparentResistivity[i] )
					/ log10 ( 1.0 + itrErr->second.error[i].first / itrData->second.Obs.apparentResistivity[i] );
				addNumDataAndSumResidual(siteID, numDataAndSumResidualEachSite, residual);
			}
			// Phase
			if( itrErr->second.error[i].second > 0.0 ){
				const double residual = ( itrData->second.Cal.phase[i] - itrData->second.Obs.phase[i] ) / itrErr->second.error[i].second;
				addNumDataAndSumResidual(siteID, numDataAndSumResidualEachSite, residual);
			}
		}
	}
	
	std::ofstream ofile( "RMS.out" );
	if( ofile.fail() )
	{
		std::cerr << "File open error : result_MT.txt !!" << std::endl;
		exit(1);
	}

	ofile << std::setw(10) << "Site" << std::setw(10) << "#Data" << std::setw(15) << "RMS" << std::endl;
	double sumResidualAll(0.0);
	int numDataAll(0);
	for( std::map<int, NumDataAndSumResidual>::const_iterator itr = numDataAndSumResidualEachSite.begin(); itr != numDataAndSumResidualEachSite.end(); ++itr ){
		const double rms = sqrt( itr->second.sumResidual / itr->second.numData );
		ofile << std::setw(10) << itr->first;
		ofile << std::setw(10) << itr->second.numData;
		ofile << std::setw(15) << std::scientific << std::setprecision(6) << rms << std::endl;
		sumResidualAll += itr->second.sumResidual;
		numDataAll += itr->second.numData;
	}
	const double rms = sqrt( sumResidualAll / numDataAll);
	ofile << std::setw(10) << "Total";
	ofile << std::setw(10) << numDataAll;
	ofile << std::setw(15) << std::scientific << std::setprecision(6) << rms << std::endl;
	ofile.close();

}

std::string convertSiteIDToSiteName( const int siteID ){

    std::map<int, std::string>::iterator itr = m_SiteIDToSiteName.find(siteID);
 
    if( itr != m_SiteIDToSiteName.end() ) {
		// found
		return itr->second;
    } else {
		// not found
		std::ostringstream oss;
		oss << siteID;
		return oss.str();
    }

}

void writeResult(){

	writeResultMT();
	writeResultVTF();
	writeResultNMT2();
	writeResultAppResAndPhs();
	writeResultNMT2AppResAndPhs();

}

void writeResultMT(){

	if( m_MTDataListAll.empty() ){
		return;
	}

	// Sort
	std::sort(m_MTDataListAll.begin(), m_MTDataListAll.end(), pairCompareMTData);

	std::string inputFileName;
	std::string delim = "";
	if( m_outputCSV ){
		inputFileName = "result_MT.csv";
		delim = ",";
	}else{
		inputFileName = "result_MT.txt";
	}
	std::ofstream ofile( inputFileName.c_str() );
	if( ofile.fail() )
	{
		std::cerr << "File open error : " << inputFileName << " !!" << std::endl;
		exit(1);
	}

	// Write header
	writeToOfstream( ofile, 10, "Site", m_outputCSV );
	writeToOfstream( ofile, 15, "Frequency", m_outputCSV );
	if( m_isImpedanceTensorConvertedToAppResAndPhase ){
		for( int i = 0; i < 4; ++i ){
			std::string app = "AppR"  + m_componentIndex[i] + "Cal";
			std::string phs = "Phs"  + m_componentIndex[i] + "Cal";
			writeToOfstream( ofile, 15, app, m_outputCSV );
			writeToOfstream( ofile, 15, phs, m_outputCSV );
		}
		if( m_typeOfDistortion != NO_DISTORTION ){
			for( int i = 0; i < 4; ++i ){
				std::string app = "AppR"  + m_componentIndex[i] + "Undist";
				std::string phs = "Phs"  + m_componentIndex[i] + "Undist";
				writeToOfstream( ofile, 15, app, m_outputCSV );
				writeToOfstream( ofile, 15, phs, m_outputCSV );
			}
		}
		for( int i = 0; i < 4; ++i ){
			std::string app = "AppR"  + m_componentIndex[i] + "Obs";
			std::string phs = "Phs"  + m_componentIndex[i] + "Obs";
			writeToOfstream( ofile, 15, app, m_outputCSV );
			writeToOfstream( ofile, 15, phs, m_outputCSV );
		}
		for( int i = 0; i < 4; ++i ){
			std::string app = "AppR"  + m_componentIndex[i] + "Err";
			std::string phs = "Phs"  + m_componentIndex[i] + "Err";
			writeToOfstream( ofile, 15, app, m_outputCSV );
			writeToOfstream( ofile, 15, phs, m_outputCSV );
		}
	}
	else{
		for( int i = 0; i < 4; ++i ){
			std::string re = "ReZ"  + m_componentIndex[i] + "Cal";
			std::string im = "ImZ"  + m_componentIndex[i] + "Cal";
			writeToOfstream( ofile, 15, re, m_outputCSV );
			writeToOfstream( ofile, 15, im, m_outputCSV );
		}
		if( m_typeOfDistortion != NO_DISTORTION ){
			for( int i = 0; i < 4; ++i ){
				std::string re = "ReZ"  + m_componentIndex[i] + "Undist";
				std::string im = "ImZ"  + m_componentIndex[i] + "Undist";
				writeToOfstream( ofile, 15, re, m_outputCSV );
				writeToOfstream( ofile, 15, im, m_outputCSV );
			}
		}
		for( int i = 0; i < 4; ++i ){
			std::string re = "ReZ"  + m_componentIndex[i] + "Obs";
			std::string im = "ImZ"  + m_componentIndex[i] + "Obs";
			writeToOfstream( ofile, 15, re, m_outputCSV );
			writeToOfstream( ofile, 15, im, m_outputCSV );
		}
		for( int i = 0; i < 4; ++i ){
			std::string re = "ReZ"  + m_componentIndex[i] + "Err";
			std::string im = "ImZ"  + m_componentIndex[i] + "Err";
			writeToOfstream( ofile, 15, re, m_outputCSV );
			writeToOfstream( ofile, 15, im, m_outputCSV );
		}
	}
	ofile << std::endl;

	ofile << std::fixed << std::scientific << std::setprecision(6);

	const double mu0 = 4.0 * M_PI * 1.0e-7;
	const double rad2deg = 180.0 / M_PI;

	double* appCal = new double[4];
	double* phsCal = new double[4];
	double* appObs = new double[4];
	double* phsObs = new double[4];
	double* appErr = new double[4];
	double* phsErr = new double[4];
	double* appCalUndist = new double[4];
	double* phsCalUndist = new double[4];
	for( std::vector<  std::pair<int, MTData> >::const_iterator itr = m_MTDataListAll.begin(); itr != m_MTDataListAll.end(); ++itr ){
		const int siteID = itr->first;
		const double freq = itr->second.freq;
		const double omega = 2.0 * M_PI * freq;
		
		// Calculate impedance tensor without galvanic distortion
		ImpedanceTensor ZUndist;
		calcUndistortedImpedanceTensor( siteID, itr->second.Cal, ZUndist );

		for( int i = 0; i < 4; ++i ){
			double error = -1;
			if( m_readTrueErrorFile ){
				std::vector< std::pair<int, MTTrueError> >::const_iterator itrErr = getIteratorToMTTrueError( siteID, freq );
				isSameError( siteID, freq, i, itrErr->second.error[i].first, itrErr->second.error[i].second );
				if( itrErr->second.error[i].first > 0.0 && itrErr->second.error[i].second > 0.0 ){
					error = std::max( itrErr->second.error[i].first, itrErr->second.error[i].second );
				}else{
					error = 1.0e10;
				}
			}else{
				isSameError( siteID, freq, i, itr->second.Err.Z[i].real(), itr->second.Err.Z[i].imag() );
				error = std::max( itr->second.Err.Z[i].real(), itr->second.Err.Z[i].imag() );
			}
			error *= FACTOR; 
			appCal[i] = std::norm(itr->second.Cal.Z[i]) / ( omega * mu0 );
			appObs[i] = std::norm(itr->second.Obs.Z[i]) / ( omega * mu0 );
			appErr[i] = 2.0 * std::abs(itr->second.Obs.Z[i]) * error / ( omega * mu0 );
			phsCal[i] = rad2deg * atan2(itr->second.Cal.Z[i].imag(), itr->second.Cal.Z[i].real() );
			phsObs[i] = rad2deg * atan2(itr->second.Obs.Z[i].imag(), itr->second.Obs.Z[i].real() );
			const double tmp = error / abs(itr->second.Obs.Z[i]);
			if( tmp <= 1.0 ){
				phsErr[i] = rad2deg * asin(tmp);
			}else{
				phsErr[i] = 180.0;
			}
			appCalUndist[i] = std::norm(ZUndist.Z[i]) / ( omega * mu0 );
			phsCalUndist[i] = rad2deg * atan2(ZUndist.Z[i].imag(), ZUndist.Z[i].real() );
		}
		writeToOfstream( ofile, 10, convertSiteIDToSiteName(itr->first), m_outputCSV );
		writeToOfstream( ofile, 15, freq, m_outputCSV );
		if( m_isImpedanceTensorConvertedToAppResAndPhase ){
			for( int i = 0; i < 4; ++i ){
				writeToOfstream( ofile, 15, appCal[i], m_outputCSV );
				writeToOfstream( ofile, 15, phsCal[i], m_outputCSV );
			}
			if( m_typeOfDistortion != NO_DISTORTION ){
				for( int i = 0; i < 4; ++i ){
					writeToOfstream( ofile, 15, appCalUndist[i], m_outputCSV );
					writeToOfstream( ofile, 15, phsCalUndist[i], m_outputCSV );
				}
			}
			for( int i = 0; i < 4; ++i ){
				writeToOfstream( ofile, 15, appObs[i], m_outputCSV );
				writeToOfstream( ofile, 15, phsObs[i], m_outputCSV );
			}
			for( int i = 0; i < 4; ++i ){
				writeToOfstream( ofile, 15, appErr[i], m_outputCSV );
				writeToOfstream( ofile, 15, phsErr[i], m_outputCSV );
			}
		}else{
			for( int i = 0; i < 4; ++i ){
				writeToOfstream( ofile, 15, itr->second.Cal.Z[i].real(), m_outputCSV );
				writeToOfstream( ofile, 15, itr->second.Cal.Z[i].imag(), m_outputCSV );
			}
			if( m_typeOfDistortion != NO_DISTORTION ){
				for( int i = 0; i < 4; ++i ){
					writeToOfstream( ofile, 15, ZUndist.Z[i].real(), m_outputCSV );
					writeToOfstream( ofile, 15, ZUndist.Z[i].imag(), m_outputCSV );
				}
			}
			for( int i = 0; i < 4; ++i ){
				writeToOfstream( ofile, 15, itr->second.Obs.Z[i].real(), m_outputCSV );
				writeToOfstream( ofile, 15, itr->second.Obs.Z[i].imag(), m_outputCSV );
			}
			for( int i = 0; i < 4; ++i ){
				writeToOfstream( ofile, 15, itr->second.Err.Z[i].real(), m_outputCSV );
				writeToOfstream( ofile, 15, itr->second.Err.Z[i].imag(), m_outputCSV );
			}
		}
		ofile << std::endl;
	}

	delete[] appCal;
	delete[] phsCal;
	delete[] appObs;
	delete[] phsObs;
	delete[] appErr;
	delete[] phsErr;
	delete[] appCalUndist;
	delete[] phsCalUndist;
	ofile.close();
}

void writeResultVTF(){

	if( m_VTFDataListAll.empty() ){
		return;
	}

	// Sort
	std::sort(m_VTFDataListAll.begin(), m_VTFDataListAll.end(), pairCompareVTFData);

	std::string inputFileName;
	std::string delim = "";
	if( m_outputCSV ){
		inputFileName = "result_VTF.csv";
		delim = ",";
	}else{
		inputFileName = "result_VTF.txt";
	}
	std::ofstream ofile( inputFileName.c_str() );
	if( ofile.fail() )
	{
		std::cerr << "File open error : " << inputFileName << " !!" << std::endl;
		exit(1);
	}

	// Write header
	writeToOfstream( ofile, 10, "Site", m_outputCSV );
	writeToOfstream( ofile, 15, "Frequency", m_outputCSV );
	for( int i = 0; i < 2; ++i ){
		std::string re = "ReT"  + m_componentIndex[i+4] + "Cal";
		std::string im = "ImT"  + m_componentIndex[i+4] + "Cal";
		writeToOfstream( ofile, 15, re, m_outputCSV );
		writeToOfstream( ofile, 15, im, m_outputCSV );
	}
	for( int i = 0; i < 2; ++i ){
		std::string re = "ReT"  + m_componentIndex[i+4] + "Obs";
		std::string im = "ImT"  + m_componentIndex[i+4] + "Obs";
		writeToOfstream( ofile, 15, re, m_outputCSV );
		writeToOfstream( ofile, 15, im, m_outputCSV );
	}
	for( int i = 0; i < 2; ++i ){
		std::string re = "ReT"  + m_componentIndex[i+4] + "Err";
		std::string im = "ImT"  + m_componentIndex[i+4] + "Err";
		writeToOfstream( ofile, 15, re, m_outputCSV );
		writeToOfstream( ofile, 15, im, m_outputCSV );
	}
	ofile << std::endl;

	ofile << std::fixed << std::scientific << std::setprecision(6);

	const double mu0 = 4.0 * M_PI * 1.0e-7;
	const double rad2deg = 180.0 / M_PI;
	
	double* reErr = new double[2];
	double* imErr = new double[2];
	for( std::vector<  std::pair<int, VTFData> >::const_iterator itr = m_VTFDataListAll.begin(); itr != m_VTFDataListAll.end(); ++itr ){
		const int siteID = itr->first;
		const double freq = itr->second.freq;

		for( int i = 0; i < 2; ++i ){
			double error = -1;
			if( m_readTrueErrorFile ){
				std::vector< std::pair<int, VTFTrueError> >::const_iterator itrErr = getIteratorToVTFTrueError( siteID, freq );
				isSameError( siteID, freq, i, itrErr->second.error[i].first, itrErr->second.error[i].second );
				if( itrErr->second.error[i].first > 0.0 && itrErr->second.error[i].second > 0.0 ){
					error = std::max( itrErr->second.error[i].first, itrErr->second.error[i].second );
				}else{
					error = 1.0e10;
				}
			}else{
				isSameError( siteID, freq, i, itr->second.Err.TZ[i].real(), itr->second.Err.TZ[i].imag() );
				error = std::max( itr->second.Err.TZ[i].real(), itr->second.Err.TZ[i].imag() );
			}
			error *= FACTOR;
			reErr[i] = error;
			imErr[i] = error;
		}
		writeToOfstream( ofile, 10, convertSiteIDToSiteName(itr->first), m_outputCSV );
		writeToOfstream( ofile, 15, freq, m_outputCSV );
		for( int i = 0; i < 2; ++i ){
			writeToOfstream( ofile, 15, itr->second.Cal.TZ[i].real(), m_outputCSV );
			writeToOfstream( ofile, 15, itr->second.Cal.TZ[i].imag(), m_outputCSV );
		}
		for( int i = 0; i < 2; ++i ){
			writeToOfstream( ofile, 15, itr->second.Obs.TZ[i].real(), m_outputCSV );
			writeToOfstream( ofile, 15, itr->second.Obs.TZ[i].imag(), m_outputCSV );
		}
		for( int i = 0; i < 2; ++i ){
			writeToOfstream( ofile, 15, reErr[i], m_outputCSV );
			writeToOfstream( ofile, 15, imErr[i], m_outputCSV );
		}
		ofile << std::endl;
	}

	delete[] reErr;
	delete[] imErr;
	ofile.close();
}

void writeResultNMT2(){

	if( m_NMT2DataListAll.empty() ){
		return;
	}

	// Sort
	std::sort(m_NMT2DataListAll.begin(), m_NMT2DataListAll.end(), pairCompareMTData);

	std::string inputFileName;
	std::string delim = "";
	if( m_outputCSV ){
		inputFileName = "result_NMT2.csv";
		delim = ",";
	}else{
		inputFileName = "result_NMT2.txt";
	}
	std::ofstream ofile( inputFileName.c_str() );
	if( ofile.fail() )
	{
		std::cerr << "File open error : " << inputFileName << " !!" << std::endl;
		exit(1);
	}

	// Write header
	writeToOfstream( ofile, 10, "Site", m_outputCSV );
	writeToOfstream( ofile, 15, "Frequency", m_outputCSV );
	if( m_isImpedanceTensorConvertedToAppResAndPhase ){
		for( int i = 0; i < 4; ++i ){
			std::string app = "AppR"  + m_componentIndex[i] + "Cal";
			std::string phs = "Phs"  + m_componentIndex[i] + "Cal";
			writeToOfstream( ofile, 15, app, m_outputCSV );
			writeToOfstream( ofile, 15, phs, m_outputCSV );
		}
		for( int i = 0; i < 4; ++i ){
			std::string app = "AppR"  + m_componentIndex[i] + "Obs";
			std::string phs = "Phs"  + m_componentIndex[i] + "Obs";
			writeToOfstream( ofile, 15, app, m_outputCSV );
			writeToOfstream( ofile, 15, phs, m_outputCSV );
		}
		for( int i = 0; i < 4; ++i ){
			std::string app = "AppR"  + m_componentIndex[i] + "Err";
			std::string phs = "Phs"  + m_componentIndex[i] + "Err";
			writeToOfstream( ofile, 15, app, m_outputCSV );
			writeToOfstream( ofile, 15, phs, m_outputCSV );
		}
	}
	else{
		for( int i = 0; i < 4; ++i ){
			std::string re = "ReZ"  + m_componentIndex[i] + "Cal";
			std::string im = "ImZ"  + m_componentIndex[i] + "Cal";
			writeToOfstream( ofile, 15, re, m_outputCSV );
			writeToOfstream( ofile, 15, im, m_outputCSV );
		}
		for( int i = 0; i < 4; ++i ){
			std::string re = "ReZ"  + m_componentIndex[i] + "Obs";
			std::string im = "ImZ"  + m_componentIndex[i] + "Obs";
			writeToOfstream( ofile, 15, re, m_outputCSV );
			writeToOfstream( ofile, 15, im, m_outputCSV );
		}
		for( int i = 0; i < 4; ++i ){
			std::string re = "ReZ"  + m_componentIndex[i] + "Err";
			std::string im = "ImZ"  + m_componentIndex[i] + "Err";
			writeToOfstream( ofile, 15, re, m_outputCSV );
			writeToOfstream( ofile, 15, im, m_outputCSV );
		}
	}
	ofile << std::endl;

	ofile << std::fixed << std::scientific << std::setprecision(6);

	double* appCal = new double[4];
	double* phsCal = new double[4];
	double* appObs = new double[4];
	double* phsObs = new double[4];
	double* appErr = new double[4];
	double* phsErr = new double[4];
	for( std::vector<  std::pair<int, MTData> >::const_iterator itr = m_NMT2DataListAll.begin(); itr != m_NMT2DataListAll.end(); ++itr ){
		const int siteID = itr->first;
		const double freq = itr->second.freq;
		const double omega = 2.0 * M_PI * freq;
		for( int i = 0; i < 4; ++i ){
			double error(-1);
			if( m_readTrueErrorFile ){
				std::vector< std::pair<int, MTTrueError> >::const_iterator itrErr = getIteratorToNMT2TrueError( siteID, freq );
				isSameError( siteID, freq, i, itrErr->second.error[i].first, itrErr->second.error[i].second );
				if( itrErr->second.error[i].first > 0.0 && itrErr->second.error[i].second > 0.0 ){
					error = std::max( itrErr->second.error[i].first, itrErr->second.error[i].second );
				}else{
					error = 1.0e10;
				}
			}else{
				isSameError( siteID, freq, i, itr->second.Err.Z[i].real(), itr->second.Err.Z[i].imag() );
				error = std::max( itr->second.Err.Z[i].real(), itr->second.Err.Z[i].imag() );
			}
			error *= FACTOR;
			appCal[i] = std::norm(itr->second.Cal.Z[i]) / ( omega * MU0 );
			appObs[i] = std::norm(itr->second.Obs.Z[i]) / ( omega * MU0 );
			appErr[i] = 2.0 * std::abs(itr->second.Obs.Z[i]) * error / ( omega * MU0 );
			phsCal[i] = RAD2DEG * atan2(itr->second.Cal.Z[i].imag(), itr->second.Cal.Z[i].real() );
			phsObs[i] = RAD2DEG * atan2(itr->second.Obs.Z[i].imag(), itr->second.Obs.Z[i].real() );
			const double tmp = error / abs(itr->second.Obs.Z[i]);
			if( tmp <= 1.0 ){
				phsErr[i] = RAD2DEG * asin(tmp);
			}else{
				phsErr[i] = 180.0;
			}
		}
		writeToOfstream( ofile, 10, convertSiteIDToSiteName(itr->first), m_outputCSV );
		writeToOfstream( ofile, 15, freq, m_outputCSV );
		if( m_isImpedanceTensorConvertedToAppResAndPhase ){
			for( int i = 0; i < 4; ++i ){
				writeToOfstream( ofile, 15, appCal[i], m_outputCSV );
				writeToOfstream( ofile, 15, phsCal[i], m_outputCSV );
			}
			for( int i = 0; i < 4; ++i ){
				writeToOfstream( ofile, 15, appObs[i], m_outputCSV );
				writeToOfstream( ofile, 15, phsObs[i], m_outputCSV );
			}
			for( int i = 0; i < 4; ++i ){
				writeToOfstream( ofile, 15, appErr[i], m_outputCSV );
				writeToOfstream( ofile, 15, phsErr[i], m_outputCSV );
			}
		}else{
			for( int i = 0; i < 4; ++i ){
				writeToOfstream( ofile, 15, itr->second.Cal.Z[i].real(), m_outputCSV );
				writeToOfstream( ofile, 15, itr->second.Cal.Z[i].imag(), m_outputCSV );
			}
			for( int i = 0; i < 4; ++i ){
				writeToOfstream( ofile, 15, itr->second.Obs.Z[i].real(), m_outputCSV );
				writeToOfstream( ofile, 15, itr->second.Obs.Z[i].imag(), m_outputCSV );
			}
			for( int i = 0; i < 4; ++i ){
				writeToOfstream( ofile, 15, itr->second.Err.Z[i].real(), m_outputCSV );
				writeToOfstream( ofile, 15, itr->second.Err.Z[i].imag(), m_outputCSV );
			}
		}
		ofile << std::endl;
	}

	delete[] appCal;
	delete[] phsCal;
	delete[] appObs;
	delete[] phsObs;
	delete[] appErr;
	delete[] phsErr;
	ofile.close();
}

void writeResultAppResAndPhs(){

	if( m_AppResAndPhsDataListAll.empty() ){
		return;
	}

	// Sort
	std::sort(m_AppResAndPhsDataListAll.begin(), m_AppResAndPhsDataListAll.end(), pairCompareAppResAndPhsData);

	std::string inputFileName;
	std::string delim = "";
	if( m_outputCSV ){
		inputFileName = "result_APP_RES_AND_PHS.csv";
		delim = ",";
	}else{
		inputFileName = "result_APP_RES_AND_PHS.txt";
	}
	std::ofstream ofile( inputFileName.c_str() );
	if( ofile.fail() )
	{
		std::cerr << "File open error : " << inputFileName << " !!" << std::endl;
		exit(1);
	}

	// Write header
	writeToOfstream( ofile, 10, "Site", m_outputCSV );
	writeToOfstream( ofile, 15, "Frequency", m_outputCSV );
	for( int i = 0; i < 4; ++i ){
		std::string app = "AppR"  + m_componentIndex[i] + "Cal";
		std::string phs = "Phs"  + m_componentIndex[i] + "Cal";
		writeToOfstream( ofile, 15, app, m_outputCSV );
		writeToOfstream( ofile, 15, phs, m_outputCSV );
	}
	if( m_typeOfDistortion != NO_DISTORTION ){
		for( int i = 0; i < 4; ++i ){
			std::string app = "AppR"  + m_componentIndex[i] + "Undist";
			std::string phs = "Phs"  + m_componentIndex[i] + "Undist";
			writeToOfstream( ofile, 15, app, m_outputCSV );
			writeToOfstream( ofile, 15, phs, m_outputCSV );
		}
	}
	for( int i = 0; i < 4; ++i ){
		std::string app = "AppR"  + m_componentIndex[i] + "Obs";
		std::string phs = "Phs"  + m_componentIndex[i] + "Obs";
		writeToOfstream( ofile, 15, app, m_outputCSV );
		writeToOfstream( ofile, 15, phs, m_outputCSV );
	}
	for( int i = 0; i < 4; ++i ){
		std::string app = "AppR"  + m_componentIndex[i] + "Err";
		std::string phs = "Phs"  + m_componentIndex[i] + "Err";
		writeToOfstream( ofile, 15, app, m_outputCSV );
		writeToOfstream( ofile, 15, phs, m_outputCSV );
	}
	ofile << std::endl;

	ofile << std::fixed << std::scientific << std::setprecision(6);
		
	double appCal[4] = { 0.0, 0.0, 0.0, 0.0 };
	double phsCal[4] = { 0.0, 0.0, 0.0, 0.0 };
	double appUndist[4] = { 0.0, 0.0, 0.0, 0.0 };
	double phsUndist[4] = { 0.0, 0.0, 0.0, 0.0 };
	double appObs[4] = { 0.0, 0.0, 0.0, 0.0 };
	double phsObs[4] = { 0.0, 0.0, 0.0, 0.0 };
	double appErr[4] = { 0.0, 0.0, 0.0, 0.0 };
	double phsErr[4] = { 0.0, 0.0, 0.0, 0.0 };
	for( std::vector<  std::pair<int, ApparentResistivityAndPhaseData> >::const_iterator itr = m_AppResAndPhsDataListAll.begin();
		itr != m_AppResAndPhsDataListAll.end(); ++itr ){
		const int siteID = itr->first;
		const double freq = itr->second.freq;
		
		// Calculate appanrent resistivity and phase without galvanic distortion
		ApparentResistivityAndPhase AppResAndPhsUndist;
		calcUndistortedApparentResistiivtyAndPhase( siteID, freq, itr->second.Cal, AppResAndPhsUndist );

		for( int i = 0; i < 4; ++i ){
			double errorAppRes(-1);
			double errorPhs(-1);
			if( m_readTrueErrorFile ){
				std::vector< std::pair<int, MTTrueError> >::const_iterator itrErr = getIteratorToAppResAndPhsTrueError( siteID, freq );
				if( itrErr->second.error[i].first > 0.0 ){
					errorAppRes =itrErr->second.error[i].first;
				}else{
					errorAppRes = 1.0e10;
				}
				if( itrErr->second.error[i].second > 0.0 ){
					errorPhs = itrErr->second.error[i].second;
				}else{
					errorPhs = 1.0e10;
				}
			}else{
				errorAppRes = itr->second.Err.apparentResistivity[i];
				errorPhs = itr->second.Err.phase[i];
			}
			appCal[i] = itr->second.Cal.apparentResistivity[i];
			appObs[i] = itr->second.Obs.apparentResistivity[i];
			appErr[i] = errorAppRes * FACTOR;
			phsCal[i] = itr->second.Cal.phase[i];
			phsObs[i] = itr->second.Obs.phase[i];
			phsErr[i] = errorPhs * FACTOR;
			if( phsErr[i] > 180.0 ){
				phsErr[i] = 180.0;
			}
			appUndist[i] = AppResAndPhsUndist.apparentResistivity[i];
			phsUndist[i] = AppResAndPhsUndist.phase[i];
		}
		writeToOfstream( ofile, 10, convertSiteIDToSiteName(itr->first), m_outputCSV );
		writeToOfstream( ofile, 15, freq, m_outputCSV );
		for( int i = 0; i < 4; ++i ){
			writeToOfstream( ofile, 15, appCal[i], m_outputCSV );
			writeToOfstream( ofile, 15, phsCal[i], m_outputCSV );
		}
		if( m_typeOfDistortion != NO_DISTORTION ){
			for( int i = 0; i < 4; ++i ){
				writeToOfstream( ofile, 15, appUndist[i], m_outputCSV );
				writeToOfstream( ofile, 15, phsUndist[i], m_outputCSV );
			}
		}
		for( int i = 0; i < 4; ++i ){
			writeToOfstream( ofile, 15, appObs[i], m_outputCSV );
			writeToOfstream( ofile, 15, phsObs[i], m_outputCSV );
		}
		for( int i = 0; i < 4; ++i ){
			writeToOfstream( ofile, 15, appErr[i], m_outputCSV );
			writeToOfstream( ofile, 15, phsErr[i], m_outputCSV );
		}
		ofile << std::endl;
	}

	ofile.close();
}

void writeResultNMT2AppResAndPhs(){

	if( m_NMT2AppResAndPhsDataListAll.empty() ){
		return;
	}

	// Sort
	std::sort(m_NMT2AppResAndPhsDataListAll.begin(), m_NMT2AppResAndPhsDataListAll.end(), pairCompareAppResAndPhsData);

	std::string inputFileName;
	std::string delim = "";
	if( m_outputCSV ){
		inputFileName = "result_NMT2_APP_RES_AND_PHS.csv";
		delim = ",";
	}else{
		inputFileName = "result_NMT2_APP_RES_AND_PHS.txt";
	}
	std::ofstream ofile( inputFileName.c_str() );
	if( ofile.fail() )
	{
		std::cerr << "File open error : " << inputFileName << " !!" << std::endl;
		exit(1);
	}

	// Write header
	writeToOfstream( ofile, 10, "Site", m_outputCSV );
	writeToOfstream( ofile, 15, "Frequency", m_outputCSV );
	for( int i = 0; i < 4; ++i ){
		std::string app = "AppR"  + m_componentIndex[i] + "Cal";
		std::string phs = "Phs"  + m_componentIndex[i] + "Cal";
		writeToOfstream( ofile, 15, app, m_outputCSV );
		writeToOfstream( ofile, 15, phs, m_outputCSV );
	}
	for( int i = 0; i < 4; ++i ){
		std::string app = "AppR"  + m_componentIndex[i] + "Obs";
		std::string phs = "Phs"  + m_componentIndex[i] + "Obs";
		writeToOfstream( ofile, 15, app, m_outputCSV );
		writeToOfstream( ofile, 15, phs, m_outputCSV );
	}
	for( int i = 0; i < 4; ++i ){
		std::string app = "AppR"  + m_componentIndex[i] + "Err";
		std::string phs = "Phs"  + m_componentIndex[i] + "Err";
		writeToOfstream( ofile, 15, app, m_outputCSV );
		writeToOfstream( ofile, 15, phs, m_outputCSV );
	}
	ofile << std::endl;

	ofile << std::fixed << std::scientific << std::setprecision(6);

	double appCal[4] = { 0.0, 0.0, 0.0, 0.0 };
	double phsCal[4] = { 0.0, 0.0, 0.0, 0.0 };
	double appObs[4] = { 0.0, 0.0, 0.0, 0.0 };
	double phsObs[4] = { 0.0, 0.0, 0.0, 0.0 };
	double appErr[4] = { 0.0, 0.0, 0.0, 0.0 };
	double phsErr[4] = { 0.0, 0.0, 0.0, 0.0 };
	for( std::vector<  std::pair<int, ApparentResistivityAndPhaseData> >::const_iterator itr = m_NMT2AppResAndPhsDataListAll.begin();
		itr != m_NMT2AppResAndPhsDataListAll.end(); ++itr ){
		const int siteID = itr->first;
		const double freq = itr->second.freq;
		for( int i = 0; i < 4; ++i ){
			double errorAppRes(-1);
			double errorPhs(-1);
			if( m_readTrueErrorFile ){
				std::vector< std::pair<int, MTTrueError> >::const_iterator itrErr = getIteratorToNMT2AppResAndPhsTrueError( siteID, freq );
				if( itrErr->second.error[i].first > 0.0 ){
					errorAppRes =itrErr->second.error[i].first;
				}else{
					errorAppRes = 1.0e10;
				}
				if( itrErr->second.error[i].second > 0.0 ){
					errorPhs = itrErr->second.error[i].second;
				}else{
					errorPhs = 1.0e10;
				}
			}else{
				errorAppRes = itr->second.Err.apparentResistivity[i];
				errorPhs = itr->second.Err.phase[i];
			}
			appCal[i] = itr->second.Cal.apparentResistivity[i];
			appObs[i] = itr->second.Obs.apparentResistivity[i];
			appErr[i] = errorAppRes * FACTOR;
			phsCal[i] = itr->second.Cal.phase[i];
			phsObs[i] = itr->second.Obs.phase[i];
			phsErr[i] = errorPhs * FACTOR;
			if( phsErr[i] > 180.0 ){
				phsErr[i] = 180.0;
			}
		}
		writeToOfstream( ofile, 10, convertSiteIDToSiteName(itr->first), m_outputCSV );
		writeToOfstream( ofile, 15, freq, m_outputCSV );
		for( int i = 0; i < 4; ++i ){
			writeToOfstream( ofile, 15, appCal[i], m_outputCSV );
			writeToOfstream( ofile, 15, phsCal[i], m_outputCSV );
		}
		for( int i = 0; i < 4; ++i ){
			writeToOfstream( ofile, 15, appObs[i], m_outputCSV );
			writeToOfstream( ofile, 15, phsObs[i], m_outputCSV );
		}
		for( int i = 0; i < 4; ++i ){
			writeToOfstream( ofile, 15, appErr[i], m_outputCSV );
			writeToOfstream( ofile, 15, phsErr[i], m_outputCSV );
		}
		ofile << std::endl;
	}

	ofile.close();
}

bool pairCompareMTData( const std::pair<int, MTData>& left, const std::pair<int, MTData>& right )
{
	if( left.first == right.first ){
		return left.second.freq < right.second.freq;
	}else{
		return left.first < right.first;
	}
}

bool pairCompareVTFData( const std::pair<int, VTFData>& left, const std::pair<int, VTFData>& right )
{
	if( left.first == right.first ){
		return left.second.freq < right.second.freq;
	}else{
		return left.first < right.first;
	}
}

bool pairCompareAppResAndPhsData( const std::pair<int, ApparentResistivityAndPhaseData>& left, const std::pair<int, ApparentResistivityAndPhaseData>& right )
{
	if( left.first == right.first ){
		return left.second.freq < right.second.freq;
	}else{
		return left.first < right.first;
	}
}

void calcUndistortedImpedanceTensor( const int siteID, const ImpedanceTensor& Z, ImpedanceTensor& ZWithoutDistortion){
	
    std::map<int, DistortionMatrix>::iterator itr = m_distortionMatrixList.find(siteID);
 
    if( itr == m_distortionMatrixList.end() ) {
		// not found
		ZWithoutDistortion = Z;
    }else{
		// found
		const double det = itr->second.C[XX] * itr->second.C[YY] - itr->second.C[XY] * itr->second.C[YX];
		const double CInvXX =   itr->second.C[YY] / det;
		const double CInvXY = - itr->second.C[XY] / det;
		const double CInvYX = - itr->second.C[YX] / det;
		const double CInvYY =   itr->second.C[XX] / det;
		ZWithoutDistortion.Z[XX] = CInvXX * Z.Z[XX] + CInvXY * Z.Z[YX];
		ZWithoutDistortion.Z[XY] = CInvXX * Z.Z[XY] + CInvXY * Z.Z[YY];
		ZWithoutDistortion.Z[YX] = CInvYX * Z.Z[XX] + CInvYY * Z.Z[YX];
		ZWithoutDistortion.Z[YY] = CInvYX * Z.Z[XY] + CInvYY * Z.Z[YY];
	}

}

void calcUndistortedApparentResistiivtyAndPhase( const int siteID, const double freq, const ApparentResistivityAndPhase& appResAndPhs, ApparentResistivityAndPhase& appResAndPhsWithoutDistortion ){
	
    std::map<int, DistortionMatrix>::iterator itr = m_distortionMatrixList.find(siteID);
	const double omega = 2 * M_PI * freq;

	if( itr == m_distortionMatrixList.end() ) {
		// not found
		appResAndPhsWithoutDistortion = appResAndPhs;
    }else{
		// found
		ImpedanceTensor Dist;
		for( int i = 0; i < 4; ++i ){
			const double absZ = sqrt( appResAndPhs.apparentResistivity[i] * MU0 * omega );
			const double phsRad = appResAndPhs.phase[i] * DEG2RAD;
			Dist.Z[i] = std::complex<double>( absZ * cos(phsRad),  absZ * sin(phsRad) );
		}
		ImpedanceTensor Undist;
		calcUndistortedImpedanceTensor( siteID, Dist, Undist );
		for( int i = 0; i < 4; ++i ){
			appResAndPhsWithoutDistortion.apparentResistivity[i] = std::norm(Undist.Z[i]) / ( omega * MU0 );
			appResAndPhsWithoutDistortion.phase[i] = atan2( Undist.Z[i].imag(), Undist.Z[i].real() ) * RAD2DEG;
		}
	}

}

int getSiteTypeFromSiteID( const int siteID ){

   std::map< int, StationType >::const_iterator itr = m_SiteIDToSiteType.find(siteID);
 
    if( itr == m_SiteIDToSiteType.end() ) {
		// not found
		std::cerr << "Site ID (" << siteID << ") is not found !!" << std::endl;
		return -1;
    }else{
		// found
		return itr->second;
	}
}

void addNumDataAndSumResidual( const int siteID, std::map<int, NumDataAndSumResidual> & numDataAndSumResidual, const double residual ){

	std::map<int, NumDataAndSumResidual>::iterator itr = numDataAndSumResidual.find(siteID);
	if( itr != numDataAndSumResidual.end() ){
		// Found
		itr->second.numData += 1;
		itr->second.sumResidual += residual * residual;
	}else{
		// Not found => insert firt data
		NumDataAndSumResidual data = { 1, residual * residual };
		numDataAndSumResidual.insert( std::make_pair( siteID, data ) );
	}

}

std::vector< std::pair<int, MTTrueError> >::const_iterator getIteratorToMTTrueError( const int siteID, const double freq ){

	std::vector< std::pair<int, MTTrueError> >::const_iterator itrErr = m_MTTrueErrorListAll.begin();
	for( ; itrErr != m_MTTrueErrorListAll.end(); ++itrErr ){
		if( siteID == itrErr->first && std::fabs(freq - itrErr->second.freq) < 1.0e-6 ){
			break;
		}
	}
	if( itrErr ==  m_MTTrueErrorListAll.end() ){
		std::cerr << "Site ID =" << siteID << " and frequency =" << freq << " is not found !!" << std::endl;
		exit(1);
	}
	return itrErr;

}

std::vector< std::pair<int, VTFTrueError> >::const_iterator getIteratorToVTFTrueError( const int siteID, const double freq ){

	std::vector< std::pair<int, VTFTrueError> >::const_iterator itrErr = m_VTFTrueErrorListAll.begin();
	for( ; itrErr != m_VTFTrueErrorListAll.end(); ++itrErr ){
		if( siteID == itrErr->first && std::fabs(freq - itrErr->second.freq) < 1.0e-6 ){
			break;
		}
	}
	if( itrErr ==  m_VTFTrueErrorListAll.end() ){
		std::cerr << "Site ID =" << siteID << " and frequency =" << freq << " is not found !!" << std::endl;
		exit(1);
	}
	return itrErr;

}

std::vector< std::pair<int, MTTrueError> >::const_iterator getIteratorToNMT2TrueError( const int siteID, const double freq ){

	std::vector< std::pair<int, MTTrueError> >::const_iterator itrErr = m_NMT2TrueErrorListAll.begin();
	for( ; itrErr != m_NMT2TrueErrorListAll.end(); ++itrErr ){
		if( siteID == itrErr->first && std::fabs(freq - itrErr->second.freq) < 1.0e-6 ){
			break;
		}
	}
	if( itrErr ==  m_NMT2TrueErrorListAll.end() ){
		std::cerr << "Site ID =" << siteID << " and frequency =" << freq << " is not found !!" << std::endl;
		exit(1);
	}
	return itrErr;

}

std::vector< std::pair<int, MTTrueError> >::const_iterator getIteratorToAppResAndPhsTrueError( const int siteID, const double freq ){

	std::vector< std::pair<int, MTTrueError> >::const_iterator itrErr = m_AppResAndPhsTrueErrorListAll.begin();
	for( ; itrErr != m_AppResAndPhsTrueErrorListAll.end(); ++itrErr ){
		if( siteID == itrErr->first && std::fabs(freq - itrErr->second.freq) < 1.0e-6 ){
			break;
		}
	}
	if( itrErr ==  m_AppResAndPhsTrueErrorListAll.end() ){
		std::cerr << "Site ID =" << siteID << " and frequency =" << freq << " is not found !!" << std::endl;
		exit(1);
	}
	return itrErr;

}

std::vector< std::pair<int, MTTrueError> >::const_iterator getIteratorToNMT2AppResAndPhsTrueError( const int siteID, const double freq ){

	std::vector< std::pair<int, MTTrueError> >::const_iterator itrErr = m_NMT2AppResAndPhsTrueErrorListAll.begin();
	for( ; itrErr != m_NMT2AppResAndPhsTrueErrorListAll.end(); ++itrErr ){
		if( siteID == itrErr->first && std::fabs(freq - itrErr->second.freq) < 1.0e-6 ){
			break;
		}
	}
	if( itrErr ==  m_NMT2AppResAndPhsTrueErrorListAll.end() ){
		std::cerr << "Site ID =" << siteID << " and frequency =" << freq << " is not found !!" << std::endl;
		exit(1);
	}
	return itrErr;

}

void isSameError( const int siteID, const double freq, const int index, const double errReal, const double errImag ){

	const double eps = 1.0e-10;

	if( abs(errReal - errImag ) > eps ){
		std::cerr << "Errors of real and imaginary part are different. ";	
		std::cerr << "Station ID: " << siteID << ", ";	
		std::cerr << "Frequency [Hz]: " << freq << ", ";		
		std::cerr << "Component index: " << index << std::endl;	
		std::cerr << "Error (real part): " << errReal << std::endl;	
		std::cerr << "Error (imaginary part): " << errImag << std::endl;	
	}
	
	return;

}

void writeToOfstream( std::ofstream& ofs, const int width, const std::string& str, const bool isCSV ){
	if( isCSV){
		ofs<< str << ",";
	}else{
		ofs << std::setw(width) << str;
	}
}

void writeToOfstream( std::ofstream& ofs, const int width, const int value, const bool isCSV ){
	if( isCSV ){
		ofs<< value << ",";
	}else{
		ofs << std::setw(width) << value;
	}
}

void writeToOfstream( std::ofstream& ofs, const int width, const double value, const bool isCSV ){
	if( isCSV ){
		ofs<< value <<  ",";
	}else{
		ofs << std::setw(width) << value;
	}
}