SELECT DISTINCT f.rfam_acc, f.type, f.description
FROM taxonomy tx
INNER JOIN rfamseq rf ON rf.ncbi_id = tx.ncbi_id
INNER JOIN full_region fr ON fr.rfamseq_acc = rf.rfamseq_acc
INNER JOIN family f ON f.rfam_acc = fr.rfam_acc
WHERE tx.tax_string LIKE 'Bacteria%' AND f.type LIKE 'Cis-reg%;';
