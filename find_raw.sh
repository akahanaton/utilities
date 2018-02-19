#!/bin/bash
#--------------------------------------------------
# for id in `cat ./ratio.na.num.all.id`
# do
#     if [[ -n $(find /storage/ProhitsStorage/ -name "*$id*") ]] ; then
#         printf "\tfind %s\n" $id 1>&2 
#         find /storage/ProhitsStorage/ -name "*$id*"
#     elif [[ -n $(find /mnt/MS/IBD_Common/ -name "*$id*") ]] ; then
#         printf "\tfind %s\n" $id  1>&2 
#         find /mnt/MS/IBD_Common/ -name "*$id*"
#     else
#         printf "\t cant find %s\n" $id  1>&2 
#     fi
#     #--------------------------------------------------
#     # find /storage/ProhitsStorage/ /mnt/MS/IBD_Common/Biopsy\ Raw\ files/ -name "*UC_CoA_Severe_382*" -exec cp -rfp {} ./raw \;
#     #-------------------------------------------------- 
# done  2>find_raw.log | grep raw | grep -v 'Lavage' |sed 's/ /\\ /g' > ratio.na.num.all.id.path
#-------------------------------------------------- 

for id in Control UC_CoA UC_CoN CD_CoA CD_CoN
do
    echo $id
    for path in `grep $id ratio.na.num.all.id.path | grep IBD_Common | sed 's/ /#/g'`
    do 
        newpath=`echo $path | sed 's/#/ /g'`
        echo  "$newpath"
        cp "$newpath" raw/$id/
    done
done
