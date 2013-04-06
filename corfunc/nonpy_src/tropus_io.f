
****************************************************************************
*                                                                          *
*      GENERAL SUBROUTINES FOR TROPUS                                      *
*      (C)Steve King, ISIS Facility, Rutherford Appleton Laboratory, UK    *
*                                                                          *
****************************************************************************

      SUBROUTINE LOAD_DATA(FNAME,DATA_CHANNELS,DATA_ARRAY)

      INTEGER MaxDim
      PARAMETER (MaxDim=4096)

      character*80 :: prompts(2)
      character*80      fname
      integer*4         data_channels
      real*4            data_array(MaxDim)

       prompts(1)='ERROR: Error opening Tropus input file: FATAL'
       prompts(2)='ERROR: Error reading Tropus input file: FATAL'

      data_channels=1
      open(12,file=fname,form='formatted',status='old',err=8001)
10       read(12,*,err=8002,end=1000)data_array(data_channels)
         data_channels=data_channels+1
         goto 10
1000  close(12)
      data_channels=data_channels-1
      goto 9999

8001  call showprompt(prompts(1))
       stop

8002  call showprompt(prompts(2))
       stop

9999  return
      end


      SUBROUTINE SAVE_DATA
     &(FNAME,DATA_CHANNELS,X_DATA_ARRAY,Y_DATA_ARRAY)

      INTEGER MaxDim
      PARAMETER (MaxDim=4096)

      character*80 :: prompts(2)
      character*80      fname
      integer*4         data_channels
C      real*4            x_data_array(MaxDim),y_data_array(MaxDim)
      real*4 :: x_data_array(MaxDim)
      real*4 :: y_data_array(MaxDim)

      prompts(1)='ERROR: Error opening Tropus output file: FATAL'
      prompts(2)='ERROR: Error writing Tropus output file: FATAL'

1000  format(f12.5,' ,',e12.5)

      open(12,file=fname,form='formatted',status='unknown',err=8001)
         do i=1,data_channels
             write(12,1000,err=8002)x_data_array(i),y_data_array(i)
         end do
      close(12)
      goto 9999

8001  call showprompt(prompts(1))
       stop

8002  close(12)
       call showprompt(prompts(2))
       stop

9999  return
      end
