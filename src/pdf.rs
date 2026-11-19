//! # IMGT/StatAssembly pdf generation
//!
//! This software analyzes BAM files to determine confidence levels on immunogenetic loci
//! (specifically Immunoglobulin\[IG\] and T-cell receptor\[TR\] loci) as well as allele confidence.
//!
//! It was created and used by the IMGT Team (<https://www.imgt.org>).
//!
//! ## License
//!
//! Available under the EUPL license.
//!
//! ## Authors
//!
//! Made by Guilhem Zeitoun and the IMGT Team.
//!
//! ## Module
//!
//! PDF summary report generation (enabled with the `pdf` feature).
//! Builds a two-page PDF listing the analyzed loci and the validated alleles.
//!

use std::collections::BTreeMap;
use std::io;
use std::io::ErrorKind::InvalidInput;
use std::path::Path;

use num_format::{Locale, ToFormattedString};

use oxidize_pdf::{Color, Document, Font, Image, Page, PageTables, TableOptions};

use crate::r#struct::{Filecrea, GeneInfos, Genehit, GenesList, LocusInfos, Positionstrand};
use crate::{AUTHOR, NAME, VERSION};

/// Generates the PDF summary report (loci table + validated alleles table) in the output directory.
pub(crate) fn generatepdf(
    infos: &[LocusInfos],
    hash: &BTreeMap<LocusInfos, Vec<Genehit>>,
    outdir: &Path,
) -> io::Result<Filecrea> {
    let locusdata: Vec<Vec<String>> = infos
        .iter()
        .map(|a| {
            vec![
                a.locusinfo.to_string(),
                a.getcontig().to_string(),
                a.getpositionstart()
                    .getobasedpos()
                    .to_formatted_string(&Locale::en),
                a.getpositionend()
                    .getobasedpos()
                    .to_formatted_string(&Locale::en),
                a.getstrand().to_string(),
                a.getstatus().to_string(),
            ]
        })
        .collect();
    let hashdata: Vec<GeneInfos> = hash
        .values()
        .flat_map(|f| f.iter().map(|b| GeneInfos::from(b.geneinfo.clone())))
        .collect();
    let elem: Vec<Vec<String>> = hashdata
        .into_iter()
        .map(|a| {
            vec![
                a.getgene().to_string(),
                a.getchromosome().to_string(),
                a.getpositionstart()
                    .getobasedpos()
                    .to_formatted_string(&Locale::en),
                a.getpositionend()
                    .getobasedpos()
                    .to_formatted_string(&Locale::en),
                a.getstrand().to_string(),
                a.getstatus().to_string(),
            ]
        })
        .collect();
    // Create a new document
    let mut doc = Document::new();
    doc.set_title(format!("IMGT/StatAssembly version {} summary", VERSION));

    // Create a page
    let mut page = Page::a4();
    let black = Color::black();
    // Add text
    page.text()
        .set_font(Font::Helvetica, 24.0)
        .set_fill_color(black)
        .at(20.0, 700.0)
        .write(&format!("{} version {}", NAME.as_str(), VERSION))
        .map_err(|d| io::Error::new(InvalidInput, d))?;

    // Add graphics
    let softwareimage = Image::from_file("images/logo_software.png")
        .map_err(|d| io::Error::new(InvalidInput, d))?;
    page.add_image("Software logo", softwareimage);
    let image =
        Image::from_file("images/logo_imgt.png").map_err(|d| io::Error::new(InvalidInput, d))?;
    page.add_image("IMGT logo", image);
    page.text()
        .set_font(Font::Helvetica, 20.0)
        .set_fill_color(black)
        .at(20.0, 600.0)
        .write("Locus found")
        .map_err(|d| io::Error::new(InvalidInput, d))?;
    let tableoptions = TableOptions::default();
    page.add_quick_table(
        locusdata,
        5.0,
        400.0,
        page.content_width() - 10.0,
        Some(tableoptions),
    )
    .map_err(|d| io::Error::new(InvalidInput, d))?;
    doc.add_page(page);
    // Create a page
    let mut page = Page::a4();
    page.text()
        .set_font(Font::Helvetica, 20.0)
        .set_fill_color(black)
        .at(2.0, 700.0)
        .write("Validated alleles")
        .map_err(|d| io::Error::new(InvalidInput, d))?;
    let tableoptions = TableOptions::default();
    page.add_quick_table(
        elem,
        5.0,
        620.0,
        page.content_width() - 10.0,
        Some(tableoptions),
    )
    .map_err(|d| io::Error::new(InvalidInput, d))?;
    // Save the document
    doc.add_page(page);
    doc.set_producer(AUTHOR);
    doc.set_author(NAME.as_str());
    /* let mut permissions = encryption::Permissions::new();
    permissions.set_copy(true).set_print_high_quality(true);
    doc.set_encryption(DocumentEncryption {
        user_password: UserPassword(String::new()),
        owner_password: OwnerPassword(String::new()),
        permissions: permissions,
        strength: oxidize_pdf::document::EncryptionStrength::Aes256,
    }); */

    let file = Filecrea::createfrompath(Path::join(outdir, "summary.pdf"));
    doc.save(file.getpath())
        .map_err(|d| io::Error::new(InvalidInput, d))?;
    Ok(file)
}
