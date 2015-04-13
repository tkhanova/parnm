#pragma once

#include <stdio.h>
#include "type.h"
#include "util.h"

void ParseArgv( int argc, char** argv, char* &matrixFile );
void PrintMatrixInCRSView( const crsMatrix* const A );
void PrintMatrixInNormalView( const crsMatrix* A );
crsMatrix* UpTriangleMatrixToFullSymmetricMatrix( const crsMatrix* const A );

