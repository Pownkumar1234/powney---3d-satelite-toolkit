const fs = require('fs');
const path = require('path');

try {
    let content = fs.readFileSync('links.json', 'utf8');
    // Strip BOM if present
    if (content.charCodeAt(0) === 0xFEFF) {
        content = content.slice(1);
    }

    let links;
    try {
        links = JSON.parse(content);
    } catch (e) {
        console.error("JSON parsing failed:");
        console.error(e.message);
        // Try to debug: print first 100 chars
        console.error("First 100 chars:", content.substring(0, 100));
        process.exit(1);
    }

    if (!Array.isArray(links)) {
        if (links && links.Path) links = [links];
        else links = [];
    }

    // Filter out the reports themselves to avoid recursion noise
    links = links.filter(l => {
        const p = l.Path || "";
        if (p.includes("found_links.txt")) return false;
        if (p.includes("links.json")) return false;
        if (p.includes("https_report.md")) return false;
        return true;
    });

    let report = "# HTTPS Links Analysis Report\n\n";
    report += `Total Occurrences Found: ${links.length}\n\n`;

    const validLinks = links.filter(l => l && l.Path);

    const byFile = {};
    validLinks.forEach(l => {
        const p = l.Path;
        if (!p) return;
        const relPath = path.relative(process.cwd(), p);
        if (!byFile[relPath]) byFile[relPath] = [];
        byFile[relPath].push(l);
    });

    report += "## Files Breakdown\n\n";
    const fileNames = Object.keys(byFile).sort();

    const dataExts = ['.json', '.md', '.txt', '.csv', '.xml', '.lock', '.map', '.cff'];
    const srcExts = ['.js', '.ts', '.tsx', '.jsx', '.html', '.css', '.scss', '.less'];

    let srcFiles = [];
    let otherFiles = [];

    fileNames.forEach(f => {
        const lower = f.toLowerCase();
        if (srcExts.some(ext => lower.endsWith(ext))) srcFiles.push(f);
        else otherFiles.push(f);
    });

    report += "### Source Code Files containing HTTPS Links\n";
    if (srcFiles.length === 0) report += "None.\n";
    srcFiles.forEach(f => {
        const absPath = path.resolve(f).replace(/\\/g, '/');
        const count = byFile[f].length;
        report += `- [${f}](file:///${absPath}) (${count} links)\n`;
    });

    report += "\n### Documentation/Data/Config Files\n";
    otherFiles.forEach(f => {
        const absPath = path.resolve(f).replace(/\\/g, '/');
        const count = byFile[f].length;
        report += `- [${f}](file:///${absPath}) (${count} links)\n`;
    });

    const domains = new Set();
    const urlRegex = /https:\/\/([a-zA-Z0-9.-]+)/g;
    validLinks.forEach(l => {
        const line = l.Line || "";
        let match;
        while ((match = urlRegex.exec(line)) !== null) {
            domains.add(match[1]);
        }
    });

    report += "\n## Unique Domains Found\n";
    Array.from(domains).sort().forEach(d => report += `- ${d}\n`);

    fs.writeFileSync('https_report.md', report);
    console.log("Report generated: https_report.md");

} catch (err) {
    console.error(err);
}
