%Preparamos el entorno de simulacion 
clear
clc
close all

%Parametros de la simulacion 
M = 16;             %Tamano de la constelacion 
k = log2(M);        %Numero de bits por simbolo
n = 10^4;           %Numero de bits a procesar     
numSamSym = 1;      %Numero de muestras por simbolo, oversampling factor

subPort = 64;       %Cantidad de sub-portadoras posibles [-32 ...... +31]
subUtil = 52;       %Cantidad de sub-portadoras utilizadas [-26 ..... 1 and +1 ..... +26]
zeroSubPort = (subPort - subUtil)/2;    %Subportadoras no utilizadas 


%Generamos una secuencia de bits aleatorios 
rng default
dataIn = randi([0 1],n,1);  %Vector columna con valores pseudoaleatorios

%................................Modulacion QAM............................
dataInMatrix = reshape(dataIn,length(dataIn)/k,k);  %Se convierten los datos de serie a paralelo 
dataSymbolsIn = bi2de(dataInMatrix);                %Se convierten los valores binarios a decimales 
dataMod = qammod(dataSymbolsIn,M,0);                %Se modula a M-QAM, vector columna con los simbolos 
dataMod = dataMod';                                 %Los simbolos QAM se guardan en un vector fila

%...............................Modulacion OFDM............................ 
[dataModP, ceros] = vec2mat(dataMod,subUtil);        %Convierte los datos de serie a paralelo, se agregan ceros por conveniencia  
rowsData = size(dataModP,1);                         %Numero de filas de los datos modulados en paralelo 

%Se asignan sub-portadoras a los datos
dataSub = [zeros(rowsData,zeroSubPort)] ;                       %Ceros para las sub-portadoras [-32 ..... -27]
dataSub = horzcat(dataSub, dataModP(:,[1:subUtil/2]));          %Simbolos asignados a [-26 ..... -1] 
dataSub = horzcat(dataSub, zeros(rowsData, 1));                 %Ceros para la sub-portadora [0] 
dataSub = horzcat(dataSub, dataModP(:,[subUtil/2+1:subUtil]));  %Simbolos asignados a [+1 ..... +26]
dataSub = horzcat(dataSub, zeros(rowsData, zeroSubPort-1));     %Ceros para las sub-portadoras [+27 ..... 31]
dataSub2 = dataSub;

%Se convierten los simbolos QAM a simbolos OFDM mediante IFFT
dataOFDM = ifft(dataSub')';     

%Se agrega un prefijo ciclico de 16 muestras que es equivalente a 0.8 useg 
cyPreSize = 16; 
cyPre = dataOFDM(:,[subPort-cyPreSize+1:subPort]);  %Se toman los datos de [+49 .... +64]        
dataOFDM = [cyPre dataOFDM];                        %Se agrega el prefijo ciclico a los simbolos OFDM 
sendData = reshape(dataOFDM.',1,size(dataOFDM,1)*size(dataOFDM,2)); %Datos paralelos a seriales

%...............................Demodulacion OFDM..........................
receiveData = vec2mat(sendData,cyPreSize+subPort);              %Se convierten los datos seriales a paralelos
receiveData = receiveData(:,cyPreSize+1:size(receiveData,2));   %Se remueve el prefijo ciclico 
receiveData = fft(receiveData')';                               %Se obtienen los simbolos QAM en paralelo
receiveData = [receiveData(:,[7:32]) receiveData(:,[34:59])];   %Se obtienen las sub-portadoras con informacion
receiveData = reshape(receiveData',1,size(receiveData,1)*size(receiveData,2)); %Simbolos QAM en serie 
receiveData = receiveData(:,1:length(receiveData)-ceros);       %Se eliminan los ceros agregados por conveniencia

%...............................Demodulacion QAM...........................
dataSymbolsOut = qamdemod(receiveData,M);
dataOutMatrix = de2bi(dataSymbolsOut,k);
dataOut = dataOutMatrix(:);  

%...............................Calculo de BER.............................
[numErrors,ber] = biterr(dataIn,dataOut);
fprintf('\nThe binary coding bit error rate = %5.2e, based on %d errors\n', ...
ber,numErrors)

