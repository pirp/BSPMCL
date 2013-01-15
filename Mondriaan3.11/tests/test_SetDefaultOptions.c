#include "Options.h"

void CheckOptions(struct opts Options) {
    if ( Options.LoadbalanceStrategy != Constant ||
         Options.LoadbalanceAdjust != AdjustNo ||
         Options.SplitStrategy != Alternate ||
         Options.Alternate_FirstDirection != FirstDirRow ||
         Options.SplitMethod != Simple ||
         Options.Partitioner != PartMondriaan ||
         Options.Metric != MetricLambda ||
         Options.DiscardFreeNets != FreeNetYes ||
         Options.SquareMatrix_DistributeVectorsEqual != EqVecNo ||
         Options.SquareMatrix_DistributeVectorsEqual_AddDummies != DumNo ||
         Options.SymmetricMatrix_UseSingleEntry != SingleEntNo ||
         Options.SymmetricMatrix_SingleEntryType != ETypeLower ||
         Options.Coarsening_NrVertices != 1 ||
         Options.Coarsening_MaxNrVtxInMatch != 2 ||
         Options.Coarsening_StopRatio != 0.25 ||
         Options.Coarsening_VtxMaxFractionOfWeight != 0.25 ||
         Options.Coarsening_MatchingStrategy != MatchRandom ||
         Options.Coarsening_InprodMatchingOrder != RandomOrder ||
         Options.Coarsening_FineSwitchLevel != 2 ||
         Options.KLFM_InitPart_NrRestarts != 10 ||
         Options.KLFM_InitPart_MaxNrLoops != 10 ||
         Options.KLFM_InitPart_MaxNrNoGainMoves != 10 ||
         Options.KLFM_Refine_MaxNrLoops != 10 ||
         Options.KLFM_Refine_MaxNrNoGainMoves != 10 ||
         Options.VectorPartition_Step3 != VecIncrease ||
         Options.VectorPartition_MaxNrLoops != 10 ||
         Options.VectorPartition_MaxNrGreedyImproves  != 10 ||
         Options.OutputFormat != OutputDMM ||
         Options.OutputMode != MultipleFiles ||
         Options.Seed != -1 ||
         Options.OrderPermutation != OrderNone) {
        printf("Error\n") ;
        exit(1);
    }
}

int main(int argc, char **argv) {
    struct opts Options;
    FILE *File;
    
    printf("Test SetDefaultOptions: ");
    
    SetDefaultOptions(&Options);
    SetOptionsFromFile(&Options, "test_SetDefaultOptions.defaults");
    CheckOptions(Options);
    
    File = fopen("test_SetDefaultOptions.check", "w");
    if (File != NULL) {
        ExportOptions(File, &Options);
        fclose(File);
        
        memset(&Options, 231, sizeof(struct opts));
        SetOptionsFromFile(&Options, "test_SetDefaultOptions.check");
        
        File = fopen("test_SetDefaultOptions.check2", "w");
        ExportOptions(File, &Options);
        fclose(File);
        
        CheckOptions(Options);
    }
    else {
        printf("Error\n");
        exit(1);
    }
 
    /* Check option values */

    printf("OK\n") ;
    exit(0);

} /* end main */
