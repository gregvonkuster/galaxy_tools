<%inherit file="/base.mako"/>
<%namespace file="/message.mako" import="render_msg" />

%if message:
    ${render_msg( message, 'done' )}
%endif

<div class="report">
    <div class="reportBody">
        <h3 align="center">
            All uploaded samples
        </h3>
        <table align="center" class="colored">
            %if len(samples) == 0:
                <tr>
                    <td colspan="2">
                        There are no uploaded samples
                    </td>
                </tr>
            %else:
                <tr class="header">
                    <td>Affy ID</td>
                    <td>Sample ID</td>
                    <td>Field call</td>
                    <td>Colony Loc</td>
                    <td>Collect Date</td>
                    <td>User Specimen ID</td>
                    <td>Registry ID</td>
                    <td>Depth</td>
                    <td>DNA Extract Method</td>
                    <td>DNA Concentration</td>
                    <td>Public After</td>
                    <td>% Miss</td>
                    <td>% Ref</td>
                    <td>% Alt</td>
                    <td>% Het</td>
                    <td>Geno ID</td>
                    <td>Pheno ID</td>
                    <td>Collect ID</td>
                    <td>Exper ID</td>
                    <td>Colony ID</td>
                    <td>Taxon ID</td>
                </tr>
                <% ctr = 0 %>
                %for sample in samples:
                    %if ctr % 2 == 1:
                        <tr class="odd_row">
                    %else:
                        <tr class="tr">
                    %endif
                        <td>${sample[0]}</td>
                        <td>${sample[1]}</td>
                        <td>${sample[2]}</td>
                        <td>${sample[3]}</td>
                        <td>${sample[4]}</td>
                        <td>${sample[5]}</td>
                        <td>${sample[6]}</td>
                        <td>${sample[7]}</td>
                        <td>${sample[8]}</td>
                        <td>${sample[9]}</td>
                        <td>${sample[10]}</td>
                        <td>${sample[11]}</td>
                        <td>${sample[12]}</td>
                        <td>${sample[13]}</td>
                        <td>${sample[14]}</td>
                        <td><a href="${h.url_for(controller='genotypes', action='of_sample', sort_id='default', order='default', genotype_id=sample[15], affy_id=sample[0])}">${sample[15]}</a></td>
                        <td><a href="${h.url_for(controller='phenotypes', action='of_sample', sort_id='default', order='default', phenotype_id=sample[16], affy_id=sample[0])}">${sample[16]}</a>${sample[16]}</a></td>
                        <td><a href="${h.url_for(controller='collectors', action='of_sample', sort_id='default', order='default', collector_id=sample[17], affy_id=sample[0])}">${sample[17]}</a>${sample[16]}</a></td>
                        <td><a href="${h.url_for(controller='experiments', action='of_sample', sort_id='default', order='default', experiment_id=sample[18], affy_id=sample[0])}">${sample[18]}</a></td>
                        <td><a href="${h.url_for(controller='colonies', action='of_sample', sort_id='default', order='default', colony_id=sample[19], affy_id=sample[0])}">${sample[19]}</a></td>
                        <td><a href="${h.url_for(controller='taxonomies', action='of_sample', sort_id='default', order='default', taxonomy_id=sample[20], affy_id=sample[0])}">${sample[20]}</a></td>
                    </tr>
                    <% ctr += 1 %>
                %endfor
            %endif
        </table>
    </div>
</div>

