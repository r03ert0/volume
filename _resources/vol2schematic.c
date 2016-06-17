#include <stdio.h>
#include <stdlib.h>
#include <math.h>

short swap_short(short a)
{
    short b;
    ((char*)&b)[0]=((char*)&a)[1];
    ((char*)&b)[1]=((char*)&a)[0];
    return b;
}
int swap_int(int a)
{
    int b;
    ((char*)&b)[0]=((char*)&a)[3];
    ((char*)&b)[1]=((char*)&a)[2];
    ((char*)&b)[2]=((char*)&a)[1];
    ((char*)&b)[3]=((char*)&a)[0];
    return b;
}
int main(int argc, char *argv[])
{
    unsigned char mine_0[]={  0x0A,0x00,0x09,0x53,0x63,0x68,0x65,0x6D,0x61,0x74,0x69,0x63};               // 'Schematic' (compound)
    unsigned char mine_h[]={  0x02,0x00,0x06,0x48,0x65,0x69,0x67,0x68,0x74,     0x00,0x02};               //        1. 'Height' (short)
    unsigned char mine_l[]={  0x02,0x00,0x06,0x4C,0x65,0x6E,0x67,0x74,0x68,     0x00,0x02};               //        2. 'Length' (short)
    unsigned char mine_w[]={  0x02,0x00,0x05,0x57,0x69,0x64,0x74,0x68,          0x00,0x02};               //        3. 'Width' (short)
    unsigned char mine_1[]={  0x09,0x00,0x08,0x45,0x6E,0x74,0x69,0x74,0x69,0x65,0x73,                     //        4. 'Entities' (list)
                              0x01,                                                                       //            tag: byte
                              0x00,0x00,0x00,0x00,                                                        //            length(int): 0
                              0x09,0x00,0x0C,0x54,0x69,0x6C,0x65,0x45,0x6E,0x74,0x69,0x74,0x69,0x65,0x73, //        5. 'TileEntities' (list)
                              0x01,                                                                       //            tag: byte
                              0x00,0x00,0x00,0x00,                                                        //            length (int): 0
                              0x09,0x00,0x09,0x54,0x69,0x6C,0x65,0x54,0x69,0x63,0x6B,0x73,                //        6. 'TileTicks' (list)
                              0x01,                                                                       //            tag: byte
                              0x00,0x00,0x00,0x00,                                                        //            length (int): 0
                              0x08,0x00,0x09,0x4D,0x61,0x74,0x65,0x72,0x69,0x61,0x6C,0x73,                //        7. 'Materials' (string)
                              0x00,0x05,0x41,0x6C,0x70,0x68,0x61};                                        //            'Alpha'
    unsigned char mine_dat[]={0x07,0x00,0x04,0x44,0x61,0x74,0x61,                                         //        8. 'Data' (byte array)
                              0x00,0x00,0x00,0x08};                                                       //            length: 8, followed by data
    unsigned char mine_bio[]={0x07,0x00,0x06,0x42,0x69,0x6F,0x6D,0x65,0x73,                               //        9. 'Biomes' (byte array)
                              0x00,0x00,0x00,0x04};                                                       //            length: 4, followed by data, 0C for example
    unsigned char mine_blk[]={0x07,0x00,0x06,0x42,0x6C,0x6F,0x63,0x6B,0x73,                               //        10. 'Blocks' (byte array)
                              0x00,0x00,0x00,0x08};                                                       //            length: 8 (=2*2*2), followed by data
    unsigned char mine_2[]=  {0x00};                                                                      //    tag End
    int h=20,l=20,w=20,nblocks;
    int i;
    int x,y,z;
    unsigned char *blocks,*sch;

    blocks=(unsigned char*)calloc(h*l*w,sizeof(char));
    nblocks=0;
    for(x=0;x<h;x++)
    for(y=0;y<l;y++)
    for(z=0;z<w;z++)
    if(pow(x-h/2,2)+pow(y-l/2,2)+pow(z-w/2,2)<8*8)
    {
        blocks[x+y*h+z*h*l]=1;
        nblocks++;
    }

    int mine_0_length=sizeof(mine_0);
    int mine_h_length=sizeof(mine_h);
    int mine_l_length=sizeof(mine_l);
    int mine_w_length=sizeof(mine_w);
    int mine_1_length=sizeof(mine_1);
    int mine_dat_length=sizeof(mine_dat);
    int mine_bio_length=sizeof(mine_bio);
    int mine_blk_length=sizeof(mine_blk);
    int mine_2_length=sizeof(mine_2);

    printf("%i %i %i\n",mine_0_length,mine_h_length,mine_l_length);
    
    int length= mine_0_length+mine_h_length+mine_l_length+mine_w_length+
                mine_1_length+mine_dat_length+mine_bio_length+mine_blk_length+mine_2_length+
                l*w+2*h*l*w;
    sch=calloc(length,sizeof(char));
    printf("%i\n",length);

    // add header (mine_0)
    length=0;
    for(i=0;i<mine_0_length;i++)
        sch[length+i]=mine_0[i];
    length+=i;

    // add height
    for(i=0;i<mine_h_length;i++)
        sch[length+i]=mine_h[i];
    ((short*)&(sch[length+9]))[0]=swap_short(h);
    length+=i;

    // add length
    for(i=0;i<mine_l_length;i++)
        sch[length+i]=mine_l[i];
    ((short*)&(sch[length+9]))[0]=swap_short(l);
    length+=i;

    // add width
    for(i=0;i<mine_w_length;i++)
        sch[length+i]=mine_w[i];
    ((short*)&(sch[length+8]))[0]=swap_short(w);
    length+=i;

    // entities (mine_1)
    for(i=0;i<mine_1_length;i++)
        sch[length++]=mine_1[i];

    // data
    for(i=0;i<mine_dat_length;i++)
        sch[length+i]=mine_dat[i];
    ((int*)&(sch[length+7]))[0]=swap_int(h*w*l);
    length+=i;

    // data (data, h*l*w 0s)
    length+=h*l*w;

    // biomes
    for(i=0;i<mine_bio_length;i++)
        sch[length+i]=mine_bio[i];
    ((int*)&(sch[length+9]))[0]=swap_int(l*w);
    length+=i;

    // biomes (data, l*w 0Cs)
    for(i=0;i<l*w;i++)
        sch[length++]=0x0C;

    // blocks
    for(i=0;i<mine_blk_length;i++)
        sch[length+i]=mine_blk[i];
    ((int*)&(sch[length+9]))[0]=swap_int(h*l*w);
    length+=i;

    // blocks (data, from blocks array)
    for(i=0;i<h*l*w;i++)
        if(blocks[i]==1)
            sch[length+i]=0x03;
    length+=i;

    // end (mine_2)
    sch[length]=mine_2[0];
    length++;

    printf("%i\n",length);

    FILE *f;
    f=fopen("myc.schematic","w");
    fwrite(sch,length,sizeof(char),f);
    fclose(f);
    
    free(blocks);
}