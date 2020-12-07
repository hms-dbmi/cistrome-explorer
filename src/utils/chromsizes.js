export const HG38_SIZES = [
    {
      "chr": "chr1",
      "size": "248956422"
    },
    {
      "chr": "chr2",
      "size": "242193529"
    },
    {
      "chr": "chr3",
      "size": "198295559"
    },
    {
      "chr": "chr4",
      "size": "190214555"
    },
    {
      "chr": "chr5",
      "size": "181538259"
    },
    {
      "chr": "chr6",
      "size": "170805979"
    },
    {
      "chr": "chr7",
      "size": "159345973"
    },
    {
      "chr": "chr8",
      "size": "145138636"
    },
    {
      "chr": "chr9",
      "size": "138394717"
    },
    {
      "chr": "chr10",
      "size": "133797422"
    },
    {
      "chr": "chr11",
      "size": "135086622"
    },
    {
      "chr": "chr12",
      "size": "133275309"
    },
    {
      "chr": "chr13",
      "size": "114364328"
    },
    {
      "chr": "chr14",
      "size": "107043718"
    },
    {
      "chr": "chr15",
      "size": "101991189"
    },
    {
      "chr": "chr16",
      "size": "90338345"
    },
    {
      "chr": "chr17",
      "size": "83257441"
    },
    {
      "chr": "chr18",
      "size": "80373285"
    },
    {
      "chr": "chr19",
      "size": "58617616"
    },
    {
      "chr": "chr20",
      "size": "64444167"
    },
    {
      "chr": "chr21",
      "size": "46709983"
    },
    {
      "chr": "chr22",
      "size": "50818468"
    },
    {
      "chr": "chrX",
      "size": "156040895"
    },
    {
      "chr": "chrY",
      "size": "57227415"
    }
];

export const HG38_START_POSITIONS = HG38_SIZES.map((d, i) => {
    return {
        chr: d.chr,
        position: HG38_SIZES.slice(0, i).reduce((acc, cur) => +acc + +cur.size, 0)
    }
});